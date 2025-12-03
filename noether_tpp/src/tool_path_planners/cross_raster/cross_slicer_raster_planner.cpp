#include <noether_tpp/tool_path_planners/cross_raster/cross_slicer_raster_planner.h>
#include <noether_tpp/utils.h>

#include <algorithm>  // std::find(), std::reverse(), std::unique()
#include <numeric>    // std::iota()
#include <stdexcept>  // std::runtime_error
#include <string>     // std::to_string()
#include <utility>    // std::move()
#include <vector>     // std::vector

#include <pcl/common/common.h>  // pcl::getMinMax3d()
#include <pcl/common/pca.h>     // pcl::PCA
#include <pcl/conversions.h>    // pcl::fromPCLPointCloud2
#include <pcl/surface/vtk_smoothing/vtk_utils.h>
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCutter.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkKdTreePointLocator.h>
#include <vtkParametricSpline.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkGenericCell.h>
#include <vtkSmartPointer.h>
#include <vtkStripper.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkVersionMacros.h>
#endif

namespace
{
static const double EPSILON = 1e-6;

Eigen::Matrix3d computeRotation(const Eigen::Vector3d& vx, const Eigen::Vector3d& vy, const Eigen::Vector3d& vz)
{
  Eigen::Matrix3d m;
  m.setIdentity();
  m.col(0) = vx.normalized();
  m.col(1) = vy.normalized();
  m.col(2) = vz.normalized();
  return m;
}

double computeLength(const vtkSmartPointer<vtkPoints>& points)
{
  const vtkIdType num_points = points->GetNumberOfPoints();
  double total_length = 0.0;
  if (num_points < 2)
  {
    return total_length;
  }

  Eigen::Vector3d p0, pf;
  for (vtkIdType i = 1; i < num_points; i++)
  {
    points->GetPoint(i - 1, p0.data());
    points->GetPoint(i, pf.data());

    total_length += (pf - p0).norm();
  }
  return total_length;
}

vtkSmartPointer<vtkPoints> enforcePointSpacing(const vtkSmartPointer<vtkPoints>& points, double total_length,
                                               double point_spacing)
{
  vtkSmartPointer<vtkPoints> new_points = vtkSmartPointer<vtkPoints>::New();

  Eigen::Vector3d a, b, current_vector;
  int i = 0;
  points->GetPoint(i, a.data());
  points->GetPoint(i + 1, b.data());
  double current_distance = (a - b).norm();
  current_vector = (b - a).normalized();
  double previous_distances = 0;

  new_points->InsertNextPoint(a.data());

  double adjusted_point_spacing = total_length / std::ceil(total_length / point_spacing);

  for (double desired_distance = point_spacing; desired_distance < total_length;
       desired_distance += adjusted_point_spacing)
  {
    while (desired_distance > previous_distances + current_distance)
    {
      previous_distances += current_distance;
      i += 1;
      points->GetPoint(i, a.data());
      points->GetPoint(i + 1, b.data());
      current_distance = (b - a).norm();
      current_vector = (b - a).normalized();
    }
    Eigen::Vector3d np = a + current_vector * (desired_distance - previous_distances);
    new_points->InsertNextPoint(np.data());
  }

  // add last point
  points->GetPoint(points->GetNumberOfPoints() - 1, b.data());
  new_points->InsertNextPoint(b.data());

  return new_points;
}

vtkIdType findClosestPoint(const Eigen::Vector3d& target, const vtkSmartPointer<vtkPoints>& points)
{
  if (points->GetNumberOfPoints() == 0)
  {
    return -1;
  }

  vtkIdType closest_id = 0;
  double min_distance = std::numeric_limits<double>::max();
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
  {
    Eigen::Vector3d pt;
    points->GetPoint(i, pt.data());
    double distance = (pt - target).norm();

    if (distance < min_distance)
    {
      min_distance = distance;
      closest_id = i;
    }
  }
  return closest_id;
}

vtkSmartPointer<vtkPoints> clipPointsAroundCenter(const vtkSmartPointer<vtkPoints>& points, vtkIdType center_idx,
                                                  double half_length)
{
  vtkSmartPointer<vtkPoints> clipped_points = vtkSmartPointer<vtkPoints>::New();

  if (center_idx < 0 || center_idx >= points->GetNumberOfPoints())
  {
    return clipped_points;
  }

  Eigen::Vector3d center;
  points->GetPoint(center_idx, center.data());

  // Store backward points in order (will be added in reverse later)
  std::vector<Eigen::Vector3d> backward_points;
  double accumulated_length = 0.0;

  // Traverse backwards from center
  for (vtkIdType i = center_idx - 1; i >= 0 && accumulated_length < half_length; i--)
  {
    Eigen::Vector3d p_curr, p_next;
    points->GetPoint(i, p_curr.data());
    points->GetPoint(i + 1, p_next.data());

    double segment_length = (p_next - p_curr).norm();

    if (accumulated_length + segment_length <= half_length)
    {
      backward_points.push_back(p_curr);
      accumulated_length += segment_length;
    }
    else
    {
      // Interpolate to get exact point at half_length
      double remaining = half_length - accumulated_length;
      Eigen::Vector3d direction = (p_curr - p_next).normalized();
      Eigen::Vector3d interpolated = p_next + direction * remaining;
      backward_points.push_back(interpolated);
      break;
    }
  }

  // Add backward points in reverse order (so they go from start to center)
  for (auto it = backward_points.rbegin(); it != backward_points.rend(); ++it)
  {
    clipped_points->InsertNextPoint(it->data());
  }

  // Add center point
  clipped_points->InsertNextPoint(center.data());

  // Traverse forward from center
  accumulated_length = 0.0;
  for (vtkIdType i = center_idx + 1; i < points->GetNumberOfPoints() && accumulated_length < half_length; i++)
  {
    Eigen::Vector3d p_curr, p_prev;
    points->GetPoint(i, p_curr.data());
    points->GetPoint(i - 1, p_prev.data());

    double segment_length = (p_curr - p_prev).norm();

    if (accumulated_length + segment_length <= half_length)
    {
      clipped_points->InsertNextPoint(p_curr.data());
      accumulated_length += segment_length;
    }
    else
    {
      // Interpolate to get exact point at half_length
      double remaining = half_length - accumulated_length;
      Eigen::Vector3d direction = (p_curr - p_prev).normalized();
      Eigen::Vector3d interpolated = p_prev + direction * remaining;
      clipped_points->InsertNextPoint(interpolated.data());
      break;
    }
  }

  return clipped_points;
}

vtkSmartPointer<vtkPoints> generateDiagonalPoints(const Eigen::Vector3d& start_corner,
                                                  const Eigen::Vector3d& end_corner, int num_points)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  if (num_points < 2)
  {
    num_points = 2;  // Minimum 2 points (start and end)
  }

  Eigen::Vector3d diagonal_vector = end_corner - start_corner;

  for (int i = 0; i < num_points; i++)
  {
    double t = static_cast<double>(i) / (num_points - 1);  // Parameter from 0 to 1
    Eigen::Vector3d point = start_corner + t * diagonal_vector;
    points->InsertNextPoint(point.data());
  }

  return points;
}

/**
 * @brief removes points that appear in multiple lists such that only one instance of that point
 *        index remains
 * @param points_lists
 */
void removeRedundant(std::vector<std::vector<vtkIdType>>& points_lists)
{
  using IdList = std::vector<vtkIdType>;
  if (points_lists.size() < 2)
  {
    return;
  }

  std::vector<std::vector<vtkIdType>> new_points_lists;
  new_points_lists.push_back(points_lists.front());
  for (std::size_t i = 1; i < points_lists.size(); i++)
  {
    IdList& current_list = points_lists[i];
    IdList new_list;
    IdList all_ids;

    // create list of all ids
    for (auto& ref_list : new_points_lists)
    {
      all_ids.insert(all_ids.end(), ref_list.begin(), ref_list.end());
    }

    for (auto& id : current_list)
    {
      // add if not found in any of the previous lists
      if (std::find(all_ids.begin(), all_ids.end(), id) == all_ids.end())
      {
        new_list.push_back(id);
      }
    }

    // add if it has enough points
    if (new_list.size() > 0)
    {
      new_points_lists.push_back(new_list);
    }
  }

  points_lists.clear();
  points_lists.assign(new_points_lists.begin(), new_points_lists.end());
}

void mergeRasterSegments(const vtkSmartPointer<vtkPoints>& points, double merge_dist,
                         std::vector<std::vector<vtkIdType>>& points_lists)
{
  using namespace Eigen;
  using IdList = std::vector<vtkIdType>;
  if (points_lists.size() < 2)
  {
    return;
  }

  std::vector<IdList> new_points_lists;
  IdList merged_list_ids;
  IdList merged_list;

  auto do_merge = [&points](const IdList& current_list, const IdList& next_list, double merge_dist,
                            IdList& merged_list) {
    Vector3d cl_point, nl_point;

    // checking front and back end points respectively
    points->GetPoint(current_list.front(), cl_point.data());
    points->GetPoint(next_list.back(), nl_point.data());
    double d = (cl_point - nl_point).norm();
    if (d < merge_dist)
    {
      merged_list.assign(next_list.begin(), next_list.end());
      merged_list.insert(merged_list.end(), current_list.begin(), current_list.end());
      return true;
    }

    // checking back and front end points respectively
    points->GetPoint(current_list.back(), cl_point.data());
    points->GetPoint(next_list.front(), nl_point.data());
    d = (cl_point - nl_point).norm();
    if (d < merge_dist)
    {
      merged_list.assign(current_list.begin(), current_list.end());
      merged_list.insert(merged_list.end(), next_list.begin(), next_list.end());
      return true;
    }
    return false;
  };

  for (std::size_t i = 0; i < points_lists.size(); i++)
  {
    if (std::find(merged_list_ids.begin(), merged_list_ids.end(), i) != merged_list_ids.end())
    {
      // already merged
      continue;
    }

    IdList current_list = points_lists[i];
    Vector3d cl_point, nl_point;
    bool seek_adjacent = true;
    while (seek_adjacent)
    {
      seek_adjacent = false;
      for (std::size_t j = i + 1; j < points_lists.size(); j++)
      {
        if (std::find(merged_list_ids.begin(), merged_list_ids.end(), j) != merged_list_ids.end())
        {
          // already merged
          continue;
        }

        merged_list.clear();
        IdList next_list = points_lists[j];
        if (do_merge(current_list, next_list, merge_dist, merged_list))
        {
          current_list = merged_list;
          merged_list_ids.push_back(static_cast<vtkIdType>(j));
          seek_adjacent = true;
          continue;
        }

        std::reverse(next_list.begin(), next_list.end());
        if (do_merge(current_list, next_list, merge_dist, merged_list))
        {
          current_list = merged_list;
          merged_list_ids.push_back(static_cast<vtkIdType>(j));
          seek_adjacent = true;
          continue;
        }
      }
    }
    new_points_lists.push_back(current_list);
  }
  points_lists.clear();
  std::copy_if(new_points_lists.begin(), new_points_lists.end(), std::back_inserter(points_lists),
               [](const IdList& l) { return l.size() > 1; });
}

noether::ToolPaths convertToPoses(const std::vector<CrossRasterConstructData>& rasters_data)
{
  noether::ToolPaths rasters_array;
  for (const CrossRasterConstructData& rd : rasters_data)
  {
    noether::ToolPath raster_path;
    std::vector<vtkSmartPointer<vtkPolyData>> raster_segments;
    raster_segments.assign(rd.raster_segments.begin(), rd.raster_segments.end());

    for (const vtkSmartPointer<vtkPolyData>& polydata : raster_segments)
    {
      noether::ToolPathSegment raster_path_segment;
      std::size_t num_points = polydata->GetNumberOfPoints();
      Eigen::Vector3d p, p_next, vx, vy, vz;
      Eigen::Isometry3d pose;
      std::vector<int> indices(num_points);
      std::iota(indices.begin(), indices.end(), 0);
      for (std::size_t i = 0; i < indices.size() - 1; i++)
      {
        int idx = indices[i];
        int idx_next = indices[i + 1];
        polydata->GetPoint(idx, p.data());
        polydata->GetPoint(idx_next, p_next.data());
        polydata->GetPointData()->GetNormals()->GetTuple(idx, vz.data());
        vx = (p_next - p).normalized();
        vy = vz.cross(vx).normalized();
        vx = vy.cross(vz).normalized();
        pose = Eigen::Translation3d(p) * Eigen::AngleAxisd(computeRotation(vx, vy, vz));
        raster_path_segment.push_back(pose);
      }

      // adding last pose
      pose.translation() = p_next;  // orientation stays the same as previous
      raster_path_segment.push_back(pose);

      raster_path.push_back(raster_path_segment);
    }

    if (!raster_path.empty())
      rasters_array.push_back(raster_path);
  }

  return rasters_array;
}

bool insertNormals(const double search_radius, vtkSmartPointer<vtkPolyData>& mesh_data_,
                   vtkSmartPointer<vtkKdTreePointLocator>& kd_tree_, vtkSmartPointer<vtkPolyData>& data,
                   vtkSmartPointer<vtkCellLocator>& cell_locator_)
{
  // Find closest cell to each point and uses its normal vector
  vtkSmartPointer<vtkDoubleArray> new_normals = vtkSmartPointer<vtkDoubleArray>::New();
  new_normals->SetNumberOfComponents(3);
  new_normals->SetNumberOfTuples(data->GetPoints()->GetNumberOfPoints());

  // get normal data
  vtkSmartPointer<vtkDataArray> normal_data = mesh_data_->GetPointData()->GetNormals();

  if (!normal_data)
  {
    return false;
  }

  Eigen::Vector3d normal_vect = Eigen::Vector3d::UnitZ();
  for (int i = 0; i < data->GetPoints()->GetNumberOfPoints(); ++i)
  {
    // locate closest cell
    Eigen::Vector3d query_point;
    data->GetPoints()->GetPoint(i, query_point.data());

    if (search_radius > 0.0)
    {
      vtkSmartPointer<vtkIdList> id_list = vtkSmartPointer<vtkIdList>::New();
      kd_tree_->FindPointsWithinRadius(search_radius, query_point.data(), id_list);
      if (id_list->GetNumberOfIds() < 1)
      {
        kd_tree_->FindClosestNPoints(1, query_point.data(), id_list);

        if (id_list->GetNumberOfIds() < 1)
        {
          return false;
        }
      }

      // compute normal average
      normal_vect = Eigen::Vector3d::Zero();
      std::size_t num_normals = 0;
      for (auto p = 0; p < id_list->GetNumberOfIds(); p++)
      {
        Eigen::Vector3d temp_normal, query_point, closest_point;
        vtkIdType p_id = id_list->GetId(p);

        if (p_id < 0)
        {
          continue;
        }

        // get normal and add it to average
        normal_data->GetTuple(p_id, temp_normal.data());
        normal_vect += temp_normal.normalized();
        num_normals++;
      }

      normal_vect /= num_normals;
    }
    else
    {
      Eigen::Vector3d closest_point;
      vtkIdType cellId;
      int subid;
      double dist2;
      cell_locator_->FindClosestPoint(query_point.data(), closest_point.data(), cellId, subid, dist2);
      mesh_data_->GetCellData()->GetNormals()->GetTuple(cellId, normal_vect.data());
    }
    normal_vect.normalize();

    // save normal
    new_normals->SetTuple3(i, normal_vect(0), normal_vect(1), normal_vect(2));
  }
  data->GetPointData()->SetNormals(new_normals);
  return true;
}

/**
 * @brief Gets the distances (normal to the cut plane) from the cut origin to the closest and furthest corners of the
 * mesh
 * @param obb_size Column-wise matrix of the object-aligned size vectors of the mesh (from PCA), relative to the mesh
 * origin
 * @param pca_centroid Centroid of the mesh determined by PCA, relative to the mesh origin
 * @param origin Origin of the raster pattern, relative to the mesh origin
 * @param dir Raster cut plane normal, relative to the mesh origin
 * @return Tuple of (min distance, max distance)
 */
static std::tuple<double, double> getDistancesToMinMaxCuts(const Eigen::Matrix3d& obb_size,
                                                           const Eigen::Vector3d& obb_centroid,
                                                           const Eigen::Vector3d& cut_origin,
                                                           const Eigen::Vector3d& cut_normal)
{
  /* Create the 8 corners of the object-aligned bounding box using vectorization
   *
   * | x0  x1  x2 | * | 1 -1  1 -1  1 -1  1 -1 | +  | ox | = | cx0  cx1  ... |
   * | y0  y1  y2 |   | 1  1 -1 -1  1  1 -1 -1 |    | oy |   | cy0  cy1  ... |
   * | z0  z1  zz |   | 1  1  1  1 -1 -1 -1 -1 |    | oz |   | cz0  cz1  ... |
   */
  Eigen::MatrixXi corner_mask(3, 8);
  // clang-format off
  corner_mask <<
      1, -1, 1, -1, 1, -1, 1, -1,
      1, 1, -1, -1, 1, 1, -1, -1,
      1, 1, 1, 1, -1, -1, -1, -1;
  // clang-format on
  Eigen::MatrixXd corners = (obb_size / 2.0) * corner_mask.cast<double>();
  corners.colwise() += obb_centroid;

  // For each corner, compute the distance from the raster origin to a plane that passes through the OBB corner and
  // whose normal is the input raster direction. Save the largest and smallest distances to the corners
  double d_max = -std::numeric_limits<double>::max();
  double d_min = std::numeric_limits<double>::max();
  for (Eigen::Index col = 0; col < corners.cols(); ++col)
  {
    Eigen::Hyperplane<double, 3> plane(cut_normal, corners.col(col));
    double d = plane.signedDistance(cut_origin);

    // Negate the signed distance because points "in front" of the cut plane have to travel in the negative plane
    // direction to get back to the plane
    d *= -1;

    if (d > d_max)
      d_max = d;
    else if (d < d_min)
      d_min = d;
  }

  return std::make_tuple(d_min, d_max);
}

void printToolPaths(const noether::ToolPaths& tool_paths)
{
  std::cout << "\n----------------------------------------" << std::endl;
  std::cout << "TOOL PATH DEBUG: number of crosses is " << tool_paths.size() << std::endl;

  for (std::size_t raster_idx = 0; raster_idx < tool_paths.size(); ++raster_idx)
  {
    const noether::ToolPath& raster_path = tool_paths[raster_idx];

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Cross #" << raster_idx << std::endl;
    std::cout << "  Number of segments: " << raster_path.size() << std::endl;

    for (std::size_t seg_idx = 0; seg_idx < raster_path.size(); ++seg_idx)
    {
      const noether::ToolPathSegment& segment = raster_path[seg_idx];
      std::cout << "    Segment #" << seg_idx << std::endl;
      std::cout << "      Number of poses/points: " << segment.size() << std::endl;

      if (segment.empty())
      {
        std::cout << "      WARNING: Empty segment!" << std::endl;
        continue;
      }

      double segment_length = 0.0;
      for (std::size_t pose_idx = 1; pose_idx < segment.size(); ++pose_idx)
      {
        Eigen::Vector3d p_prev = segment[pose_idx - 1].translation();
        Eigen::Vector3d p_curr = segment[pose_idx].translation();
        segment_length += (p_curr - p_prev).norm();
      }
      std::cout << "      Segment length: " << segment_length << std::endl;

      // Print first pose
      Eigen::Vector3d poseInit = segment.front().translation();
      Eigen::Vector3d normalInit = segment.front().rotation().col(2);  // Z-axis is normal

      std::cout << "      First pose/point:" << std::endl;
      std::cout << "        Position: [" << poseInit.x() << ", " << poseInit.y() << ", " << poseInit.z() << "]"
                << std::endl;
      std::cout << "        Normal: [" << normalInit.x() << ", " << normalInit.y() << ", " << normalInit.z() << "]"
                << std::endl;

      // Print last pose
      Eigen::Vector3d poseEnd = segment.back().translation();
      Eigen::Vector3d normalEnd = segment.back().rotation().col(2);

      std::cout << "      Last pose:" << std::endl;
      std::cout << "        Position: [" << poseEnd.x() << ", " << poseEnd.y() << ", " << poseEnd.z() << "]"
                << std::endl;
      std::cout << "        Normal: [" << normalEnd.x() << ", " << normalEnd.y() << ", " << normalEnd.z() << "]"
                << std::endl;

      Eigen::Vector3d direction = (poseEnd - poseInit).normalized();
      std::cout << "      Direction: [" << direction.x() << ", " << direction.y() << ", " << direction.z() << "]"
                << std::endl;

      // Print all points
      std::cout << "      All points of segment #" << seg_idx << ":" << std::endl;
      for (std::size_t pose_idx = 0; pose_idx < segment.size(); ++pose_idx)
      {
        Eigen::Vector3d pos = segment[pose_idx].translation();
        std::cout << "        Pose #" << pose_idx << ": [" << pos.x() << ", " << pos.y() << ", " << pos.z() << "]"
                  << std::endl;
      }
    }
  }
}

void printCrossRasters(const vtkIdType& pt_idx, const Eigen::Vector3d& current_loc, const double& diagonal1_length,
                       const double& diagonal2_length, const std::vector<CrossRasterConstructData>& point_diagonals,
                       Eigen::Vector3d& mesh_normal)
{
  std::cout << "\n----------------------------------------" << std::endl;
  std::cout << "CROSS DEBUG: at point " << pt_idx << std::endl;
  std::cout << "  Center location: " << current_loc.transpose() << std::endl;
  std::cout << "  Diagonal 1 length: " << diagonal1_length << std::endl;
  std::cout << "  Diagonal 2 length: " << diagonal2_length << std::endl;
  // Print Diagonal 1
  if (point_diagonals.size() >= 1 && !point_diagonals[0].raster_segments.empty())
  {
    vtkSmartPointer<vtkPolyData> diag1 = point_diagonals[0].raster_segments[0];
    std::cout << "\n  DIAGONAL 1: " << diag1->GetNumberOfPoints() << " points" << std::endl;

    for (vtkIdType i = 0; i < diag1->GetNumberOfPoints(); i++)
    {
      Eigen::Vector3d pos, normal;
      diag1->GetPoint(i, pos.data());
      diag1->GetPointData()->GetNormals()->GetTuple(i, normal.data());

      std::cout << "    Point " << i << ": pos=[" << pos.transpose() << "]";
      std::cout << " normal=[" << normal.transpose() << "]";

      // Check if normal is flipped (dot product with mesh normal)
      double dot = normal.dot(mesh_normal);
      if (dot < 0)
      {
        std::cout << " *** WARNING: Normal flipped! (dot=" << dot << ") ***";
      }
      std::cout << std::endl;
    }

    // Check direction consistency
    if (diag1->GetNumberOfPoints() >= 2)
    {
      Eigen::Vector3d start, end;
      diag1->GetPoint(0, start.data());
      diag1->GetPoint(diag1->GetNumberOfPoints() - 1, end.data());
      Eigen::Vector3d direction = (end - start).normalized();

      std::cout << "    Direction vector: [" << direction.transpose() << "]" << std::endl;
      std::cout << "    Start to end distance: " << (end - start).norm() << std::endl;
    }
  }
  // Print Diagonal 2
  if (point_diagonals.size() >= 2 && !point_diagonals[1].raster_segments.empty())
  {
    vtkSmartPointer<vtkPolyData> diag2 = point_diagonals[1].raster_segments[0];
    std::cout << "\n  DIAGONAL 2: " << diag2->GetNumberOfPoints() << " points" << std::endl;

    for (vtkIdType i = 0; i < diag2->GetNumberOfPoints(); i++)
    {
      Eigen::Vector3d pos, normal;
      diag2->GetPoint(i, pos.data());
      diag2->GetPointData()->GetNormals()->GetTuple(i, normal.data());

      std::cout << "    Point " << i << ": pos=[" << pos.transpose() << "]";
      std::cout << " normal=[" << normal.transpose() << "]";

      // Check if normal is flipped (dot product with mesh normal)
      double dot = normal.dot(-mesh_normal);
      if (dot < 0)
      {
        std::cout << " WARNING: Normal flipped! (dot=" << dot << ") ***";
      }
      std::cout << std::endl;
    }

    // Check direction consistency
    if (diag2->GetNumberOfPoints() >= 2)
    {
      Eigen::Vector3d start, end;
      diag2->GetPoint(0, start.data());
      diag2->GetPoint(diag2->GetNumberOfPoints() - 1, end.data());
      Eigen::Vector3d direction = (end - start).normalized();
      std::cout << "    Direction vector: [" << direction.transpose() << "]" << std::endl;
      std::cout << "    Start to end distance: " << (end - start).norm() << std::endl;
    }
  }
}

void determinePrincipalAxes(const pcl::PolygonMesh& mesh, Eigen::Vector3d& mesh_normal, Eigen::Matrix3d& pca_vecs,
                            Eigen::Vector3d& centroid)
{
  // Use Original PCL point cloud
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

  // Perform PCA analysis
  pcl::PCA<pcl::PointXYZ> pca;
  pca.setInputCloud(cloud);

  // Get the extents of the point cloud relative to its mean and principal axes
  pcl::PointCloud<pcl::PointXYZ> proj;
  pca.project(*cloud, proj);
  pcl::PointXYZ min, max;
  pcl::getMinMax3D(proj, min, max);
  Eigen::Array3f scales = max.getArray3fMap() - min.getArray3fMap();

  mesh_normal = pca.getEigenVectors().col(2).cast<double>().normalized();
  centroid = pca.getMean().head<3>().cast<double>();
  pca_vecs = (pca.getEigenVectors().array().rowwise() * scales.transpose()).cast<double>();
  // Ensure consistency of direction of eigenvectors
  if (pca_vecs.col(0).sum() < 0)
  {
    pca_vecs.col(0) = -pca_vecs.col(0);
  }
  if (pca_vecs.col(1).sum() < 0)
  {
    pca_vecs.col(1) = -pca_vecs.col(1);
  }
  pca_vecs.col(2) = pca_vecs.col(0).cross(pca_vecs.col(1)).normalized() * scales(2);
}

void generateIntersectionData(const Eigen::Vector3d& cut_origin, const Eigen::Vector3d& cut_normal,
                              vtkSmartPointer<vtkAppendPolyData>& raster_data, vtkSmartPointer<vtkPolyData>& mesh_data)
{
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
  vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();

  plane->SetOrigin(cut_origin.x(), cut_origin.y(), cut_origin.z());
  plane->SetNormal(cut_normal.x(), cut_normal.y(), cut_normal.z());

  cutter->SetCutFunction(plane);
  cutter->SetInputData(mesh_data);
  cutter->SetSortBy(1);
  cutter->SetGenerateTriangles(false);
  cutter->Update();

  stripper->SetInputConnection(cutter->GetOutputPort());
  stripper->JoinContiguousSegmentsOn();
  stripper->SetMaximumLength(mesh_data->GetNumberOfPoints());
  stripper->Update();

  if (stripper->GetErrorCode() != vtkErrorCode::NoError)
  {
    return;
  }

  for (int r = 0; r < stripper->GetNumberOfOutputPorts(); r++)
  {
    raster_data->AddInputData(stripper->GetOutput(r));
  }
}

static std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> calculateDiagonalCuts(
    const Eigen::Vector3d& cross_center, const Eigen::Vector3d& cut_direction, const Eigen::Vector3d& cut_normal,
    const Eigen::Vector3d& mesh_normal, const double cross_width, const double cross_height)
{
  // Width direction (BL to BR)
  Eigen::Vector3d u_axis = cut_direction.normalized();
  // Height direction (BL to TR)
  Eigen::Vector3d v_axis = cut_normal.normalized();

  // Bottom-Left (BL)
  Eigen::Vector3d corner_BL = cross_center - cross_width / 2.0 * u_axis - cross_height / 2.0 * v_axis;
  // Bottom-Right (BR)
  Eigen::Vector3d corner_BR = corner_BL + cross_width * u_axis;
  // Top-Left (TL)
  Eigen::Vector3d corner_TL = corner_BL + cross_height * v_axis;
  // Top-Right (TR)
  Eigen::Vector3d corner_TR = corner_BL + cross_width * u_axis + cross_height * v_axis;

  // Diagonal 1: BL to TR
  Eigen::Vector3d cut_origin_d1 = corner_BL;
  Eigen::Vector3d cut_direction_d1 = (corner_TR - corner_BL).normalized();

  // Diagonal 2: TL to BR
  Eigen::Vector3d cut_origin_d2 = corner_TL;
  Eigen::Vector3d cut_direction_d2 = (corner_BR - corner_TL).normalized();

  return std::make_tuple(cut_origin_d1, cut_direction_d1, cut_origin_d2, cut_direction_d2);
}

}  // namespace

namespace noether
{
CrossSlicerRasterPlanner::CrossSlicerRasterPlanner(CrossDirectionGenerator::ConstPtr dir_gen,
                                                   CrossOriginGenerator::ConstPtr origin_gen)
  : CrossRasterPlanner(std::move(dir_gen), std::move(origin_gen))
{
}

void CrossSlicerRasterPlanner::setMinSegmentSize(const double min_segment_size)
{
  min_segment_size_ = min_segment_size;
}

void CrossSlicerRasterPlanner::setSearchRadius(const double search_radius)
{
  search_radius_ = search_radius;
}

void CrossSlicerRasterPlanner::generateRastersBidirectionally(const bool bidirectional)
{
  bidirectional_ = bidirectional;
}

void CrossSlicerRasterPlanner::setCrossDimensions(const double width, const double height)
{
  cross_width_ = width;
  cross_height_ = height;
  cross_length_ = std::sqrt(std::pow(cross_width_, 2) + std::pow(cross_height_, 2));
}

void CrossSlicerRasterPlanner::setCrossSpacing(const double cross_spacing)
{
  cross_spacing_ = cross_spacing;
}

ToolPaths CrossSlicerRasterPlanner::planImpl(const pcl::PolygonMesh& mesh) const
{
  vtkSmartPointer<vtkPolyData> mesh_data = updateMesh(mesh);

  // Build cell locator and kd_tree to recover normals
  vtkSmartPointer<vtkKdTreePointLocator> kd_tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
  kd_tree->SetDataSet(mesh_data);
  kd_tree->BuildLocator();

  vtkSmartPointer<vtkCellLocator> cell_locator = vtkSmartPointer<vtkCellLocator>::New();
  cell_locator->SetDataSet(mesh_data);
  cell_locator->BuildLocator();

  Eigen::Vector3d mesh_normal, cut_direction, cut_normal, cut_origin, centroid;
  Eigen::Matrix3d pca_vecs;
  computeCuttingPlaneParameters(mesh, mesh_normal, pca_vecs, centroid, cut_direction, cut_normal, cut_origin);

  auto [d_min_cut, d_max_cut] = calculateCuttingRange(pca_vecs, centroid, cut_origin, cut_normal);

  const double cut_span = std::abs(d_max_cut - d_min_cut);
  const auto num_planes = static_cast<std::size_t>(std::ceil(cut_span / line_spacing_));
  const Eigen::Vector3d start_loc = cut_origin + cut_normal * d_min_cut;

  // Generate primary rasters
  vtkSmartPointer<vtkAppendPolyData> raster_data = vtkSmartPointer<vtkAppendPolyData>::New();
  for (std::size_t i = 0; i < num_planes + 1; i++)
  {
    Eigen::Vector3d current_loc = start_loc + i * line_spacing_ * cut_normal;
    generateIntersectionData(current_loc, cut_normal, raster_data, mesh_data);
  }

  // Process raster slices into segments
  std::vector<CrossRasterConstructData> primary_rasters;
  raster_data->Update();
  vtkIdType num_slices = raster_data->GetTotalNumberOfInputConnections();
  for (std::size_t i = 0; i < num_slices; i++)
  {
    CrossRasterConstructData raster =
        processRasterSlice(raster_data->GetInput(i), mesh_data, cut_direction, kd_tree, cell_locator, i);

    // Save raster
    if (!raster.raster_segments.empty())
      primary_rasters.push_back(raster);
  }

  // Generate crosses
  std::vector<CrossRasterConstructData> cross_rasters;

  // Iterate through all primary rasters
  for (std::size_t primary_raster_idx = 0; primary_raster_idx < primary_rasters.size(); primary_raster_idx++)
  {
    // Iterate through all segments of the primary raster
    for (std::size_t i = 0; i < primary_rasters[primary_raster_idx].raster_segments.size(); i++)
    {
      vtkSmartPointer<vtkPolyData> segment = primary_rasters[primary_raster_idx].raster_segments[i];

      // Iterate through all points of the segment: each point is the center of a cross
      for (vtkIdType pt_idx = 0; pt_idx < segment->GetNumberOfPoints(); pt_idx++)
      {
        double p[3];
        segment->GetPoint(pt_idx, p);
        Eigen::Vector3d current_center(p[0], p[1], p[2]);

        std::vector<CrossRasterConstructData> diagonal_rasters = processDiagonalsAtPoint(
            current_center, mesh_data, cut_direction, cut_normal, mesh_normal, kd_tree, cell_locator);

        // Validate and add diagonal pair
        double diagonal1_length, diagonal2_length;
        if (validateDiagonalPair(diagonal_rasters, diagonal1_length, diagonal2_length))
        {
          cross_rasters.insert(cross_rasters.end(), diagonal_rasters.begin(), diagonal_rasters.end());
          printCrossRasters(pt_idx, current_center, diagonal1_length, diagonal2_length, diagonal_rasters, mesh_normal);
        }
      }
    }
  }

  // Convert to tool paths
  ToolPaths tool_paths = convertToPoses(cross_rasters);
  // printToolPaths(tool_paths);
  return tool_paths;
}

vtkSmartPointer<vtkPolyData> CrossSlicerRasterPlanner::updateMesh(const pcl::PolygonMesh& mesh) const
{
  if (!hasNormals(mesh))
  {
    std::stringstream ss;
    ss << "The input mesh does not have vertex normals, which are required for the plane slice raster tool path "
          "planner. "
       << "Use a MeshModifier to generate vertex normals for the mesh, or provide a different mesh with vertex "
          "normals.";
    throw std::runtime_error(ss.str());
  }

  // Convert input mesh to VTK type & calculate normals if necessary
  vtkSmartPointer<vtkPolyData> mesh_data = vtkSmartPointer<vtkPolyData>::New();
  pcl::VTKUtils::mesh2vtk(mesh, mesh_data);
  mesh_data->BuildLinks();
  mesh_data->BuildCells();
  if (!mesh_data->GetPointData()->GetNormals() || !mesh_data->GetCellData()->GetNormals())
  {
    vtkSmartPointer<vtkPolyDataNormals> normal_generator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_generator->SetInputData(mesh_data);
    normal_generator->ComputePointNormalsOn();
    normal_generator->SetComputeCellNormals(!mesh_data->GetCellData()->GetNormals());
    normal_generator->SetFeatureAngle(M_PI_2);
    normal_generator->SetSplitting(true);
    normal_generator->SetConsistency(true);
    normal_generator->SetAutoOrientNormals(false);
    normal_generator->SetFlipNormals(false);
    normal_generator->SetNonManifoldTraversal(false);
    normal_generator->Update();

    if (!mesh_data->GetPointData()->GetNormals())
    {
      mesh_data->GetPointData()->SetNormals(normal_generator->GetOutput()->GetPointData()->GetNormals());
    }
    if (!mesh_data->GetCellData()->GetNormals())
    {
      mesh_data->GetCellData()->SetNormals(normal_generator->GetOutput()->GetCellData()->GetNormals());
    }
  }

  return mesh_data;
}

void CrossSlicerRasterPlanner::computeCuttingPlaneParameters(const pcl::PolygonMesh& mesh, Eigen::Vector3d& mesh_normal,
                                                             Eigen::Matrix3d& pca_vecs, Eigen::Vector3d& centroid,
                                                             Eigen::Vector3d& cut_direction,
                                                             Eigen::Vector3d& cut_normal,
                                                             Eigen::Vector3d& cut_origin) const
{
  // Use principal component analysis (PCA) to determine the principal axes of the mesh
  determinePrincipalAxes(mesh, mesh_normal, pca_vecs, centroid);

  // Get the cutting plane parameters
  cut_direction = dir_gen_->generate(mesh);
  cut_normal = (cut_direction.normalized().cross(mesh_normal)).normalized();
  cut_origin = origin_gen_->generate(mesh);

  std::cout << "cut_normal: " << cut_normal << std::endl;
  std::cout << "mesh_normal: " << mesh_normal << std::endl;
  std::cout << "cut_direction: " << cut_direction << std::endl;
}

std::pair<double, double> CrossSlicerRasterPlanner::calculateCuttingRange(const Eigen::Matrix3d& pca_vecs,
                                                                          const Eigen::Vector3d& centroid,
                                                                          const Eigen::Vector3d& cut_origin,
                                                                          const Eigen::Vector3d& cut_normal) const
{
  double d_min_cut, d_max_cut;
  std::tie(d_min_cut, d_max_cut) = getDistancesToMinMaxCuts(pca_vecs, centroid, cut_origin, cut_normal);

  // Adjust range for unidirectional rasters
  // If we don't want to generate rasters bidirectionally...
  if (!bidirectional_)
  {
    // ... and the furthest cut distance is behind the cutting plane, then don't make any cuts
    if (d_max_cut < 0.0)
      return { 0.0, 0.0 };

    // ... and the closest cut distance is behined the cutting plane, set the distance to the closest cutting plane
    // equal to zero (i.e., at the nominal cut origin)
    if (d_min_cut < 0.0)
      d_min_cut = 0.0;
  }

  return { d_min_cut, d_max_cut };
}

CrossRasterConstructData CrossSlicerRasterPlanner::processRasterSlice(vtkSmartPointer<vtkPolyData> raster_lines,
                                                                      vtkSmartPointer<vtkPolyData> mesh_data,
                                                                      const Eigen::Vector3d& cut_direction,
                                                                      vtkSmartPointer<vtkKdTreePointLocator> kd_tree,
                                                                      vtkSmartPointer<vtkCellLocator> cell_locator,
                                                                      std::size_t slice_index) const
{
  CrossRasterConstructData raster;

  if (raster_lines->GetNumberOfLines() == 0)
    return raster;

  // Extract and process segment IDs
  std::vector<std::vector<vtkIdType>> raster_ids = extractRasterSegmentIds(raster_lines);

  if (raster_ids.empty())
    return raster;

  // Remove redundant indices and merge segments
  removeRedundant(raster_ids);
  mergeRasterSegments(raster_lines->GetPoints(), min_hole_size_, raster_ids);

  // Process each segment
  for (auto& point_ids : raster_ids)
  {
    vtkSmartPointer<vtkPoints> points = processAndAlignPoints(point_ids, raster_lines, cut_direction);

    // if (!points)
    //   continue;

    // compute length and add points if segment length is greater than threshold
    double line_length = ::computeLength(points);
    if (line_length > min_segment_size_ && points->GetNumberOfPoints() > 1)
    {
      // Enforce point spacing
      vtkSmartPointer<vtkPoints> new_points = enforcePointSpacing(points, line_length, cross_spacing_);

      // add new points to segment
      vtkSmartPointer<vtkPolyData> segment_data = vtkSmartPointer<vtkPolyData>::New();
      segment_data->SetPoints(new_points);

      // inserting normals
      if (!insertNormals(search_radius_, mesh_data, kd_tree, segment_data, cell_locator))
      {
        throw std::runtime_error("Could not insert normals for segment " +
                                 std::to_string(raster.raster_segments.size()) + " of raster " +
                                 std::to_string(slice_index));
      }

      // saving into raster
      raster.raster_segments.push_back(segment_data);
      raster.segment_lengths.push_back(line_length);
    }
  }

  return raster;
}

std::vector<std::vector<vtkIdType>>
CrossSlicerRasterPlanner::extractRasterSegmentIds(vtkSmartPointer<vtkPolyData> raster_lines) const
{
  std::vector<std::vector<vtkIdType>> raster_ids;

#if VTK_MAJOR_VERSION > 7
  const vtkIdType* indices;
#else
  vtkIdType* indices;
#endif
  vtkIdType num_points;
  vtkCellArray* cells = raster_lines->GetLines();

  for (cells->InitTraversal(); cells->GetNextCell(num_points, indices);)
  {
    std::vector<vtkIdType> point_ids;

    for (vtkIdType i = 0; i < num_points; i++)
    {
      if (std::find(point_ids.begin(), point_ids.end(), indices[i]) == point_ids.end())
      {
        point_ids.push_back(indices[i]);
      }
    }

    if (point_ids.empty())
      continue;

    // Remove duplicates
    auto iter = std::unique(point_ids.begin(), point_ids.end());
    point_ids.erase(iter, point_ids.end());

    raster_ids.push_back(point_ids);
  }

  return raster_ids;
}

vtkSmartPointer<vtkPoints> CrossSlicerRasterPlanner::processAndAlignPoints(const std::vector<vtkIdType>& point_ids,
                                                                           vtkSmartPointer<vtkPolyData> raster_lines,
                                                                           const Eigen::Vector3d& cut_direction) const
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // Populating with points
  std::for_each(point_ids.begin(), point_ids.end(), [&points, &raster_lines](const vtkIdType& id) {
    std::array<double, 3> p;
    raster_lines->GetPoint(id, p.data());
    points->InsertNextPoint(p.data());
  });

  // If points are not aligned with cut_direction, reverse them
  Eigen::Vector3d p0, p1;
  points->GetPoint(0, p0.data());
  points->GetPoint(points->GetNumberOfPoints() - 1, p1.data());
  if ((p1 - p0).dot(cut_direction) < 0)
  {
    vtkSmartPointer<vtkPoints> reversed_points = vtkSmartPointer<vtkPoints>::New();
    for (vtkIdType pi = points->GetNumberOfPoints() - 1; pi >= 0; pi--)
    {
      std::array<double, 3> p;
      points->GetPoint(pi, p.data());
      reversed_points->InsertNextPoint(p.data());
    }
    points = reversed_points;
    // return reversed_points;
  }

  return points;
}

std::vector<CrossRasterConstructData> CrossSlicerRasterPlanner::processDiagonalsAtPoint(
    const Eigen::Vector3d& current_loc, vtkSmartPointer<vtkPolyData> mesh_data, const Eigen::Vector3d& cut_direction,
    const Eigen::Vector3d& cut_normal, const Eigen::Vector3d& mesh_normal,
    vtkSmartPointer<vtkKdTreePointLocator> kd_tree, vtkSmartPointer<vtkCellLocator> cell_locator) const
{
  std::vector<CrossRasterConstructData> temp_diagonals;

  // Calculate diagonal cut parameters
  Eigen::Vector3d cut_origin_d1, cut_direction_d1, cut_origin_d2, cut_direction_d2;
  std::tie(cut_origin_d1, cut_direction_d1, cut_origin_d2, cut_direction_d2) =
      calculateDiagonalCuts(current_loc, cut_direction, cut_normal, mesh_normal, cross_width_, cross_height_);

  // Iterate over both cross diagonals
  for (int diagonal_idx = 0; diagonal_idx < 2; diagonal_idx++)
  {
    Eigen::Vector3d current_cut_direction = (diagonal_idx == 0) ? cut_direction_d1 : cut_direction_d2;
    Eigen::Vector3d current_cut_normal = current_cut_direction.cross(mesh_normal).normalized();

    // Generate intersection points for the diagonal
    vtkSmartPointer<vtkAppendPolyData> single_cut_data = vtkSmartPointer<vtkAppendPolyData>::New();
    generateIntersectionData(current_loc, current_cut_normal, single_cut_data, mesh_data);
    single_cut_data->Update();

    vtkIdType num_cuts = single_cut_data->GetTotalNumberOfInputConnections();
    if (num_cuts == 0)
    {
      std::cout << "  No intersections found for this cut" << std::endl;
      continue;
    }

    // Process intersection points (similar to processRasterSlice but adapted for diagonals)
    for (vtkIdType cut_idx = 0; cut_idx < num_cuts; cut_idx++)
    {
      vtkSmartPointer<vtkPolyData> raster_lines_diagonal = single_cut_data->GetInput(cut_idx);
      std::vector<std::vector<vtkIdType>> raster_ids_diagonal = extractRasterSegmentIds(raster_lines_diagonal);

      if (raster_ids_diagonal.empty())
        continue;

      removeRedundant(raster_ids_diagonal);
      mergeRasterSegments(raster_lines_diagonal->GetPoints(), min_hole_size_, raster_ids_diagonal);

      for (auto& point_ids : raster_ids_diagonal)
      {
        vtkSmartPointer<vtkPoints> points =
            processAndAlignPoints(point_ids, raster_lines_diagonal, current_cut_direction);

        if (!points)
          continue;

        // Find closest point to current_loc
        vtkIdType closest_point_id = findClosestPoint(current_loc, points);
        if (closest_point_id < 0)
          continue;

        double line_length = ::computeLength(points);
        if (line_length > min_segment_size_ && points->GetNumberOfPoints() > 1)
        {
          // Enforce point spacing
          vtkSmartPointer<vtkPoints> new_points = enforcePointSpacing(points, line_length, point_spacing_);

          // Clip points to diagonal length
          double diagonal_half_length = cross_length_ / 2.0;
          vtkIdType respaced_center_idx = findClosestPoint(current_loc, new_points);

          vtkSmartPointer<vtkPoints> clipped_points =
              clipPointsAroundCenter(new_points, respaced_center_idx, diagonal_half_length);

          if (clipped_points->GetNumberOfPoints() < 2)
            continue;

          double clipped_length = ::computeLength(clipped_points);
          if (clipped_length < min_segment_size_)
            continue;

          vtkSmartPointer<vtkPoints> final_points = enforcePointSpacing(clipped_points, clipped_length, point_spacing_);

          // Create segment with normals
          vtkSmartPointer<vtkPolyData> segment_data = vtkSmartPointer<vtkPolyData>::New();
          segment_data->SetPoints(final_points);

          if (!insertNormals(search_radius_, mesh_data, kd_tree, segment_data, cell_locator))
            continue;

          CrossRasterConstructData diagonal_raster;
          diagonal_raster.raster_segments.push_back(segment_data);
          diagonal_raster.segment_lengths.push_back(clipped_length);
          temp_diagonals.push_back(diagonal_raster);
        }
      }
    }
  }

  return temp_diagonals;
}

bool CrossSlicerRasterPlanner::validateDiagonalPair(const std::vector<CrossRasterConstructData>& diagonals,
                                                    double& diagonal1_length, double& diagonal2_length) const
{
  if (diagonals.size() != 2)
  {
    return false;
  }

  diagonal1_length = diagonals[0].segment_lengths[0];
  diagonal2_length = diagonals[1].segment_lengths[0];

  // Check minimum length requirement (90% of expected)
  const double min_length_tolerance = 0.10;
  double min_required_length = cross_length_ * (1.0 - min_length_tolerance);

  if (diagonal1_length < min_required_length || diagonal2_length < min_required_length)
  {
    // std::cout << "  Cross REJECTED: diagonal(s) too short. D1=" << diagonal1_length << ", D2=" << diagonal2_length
    //           << ", Min required=" << min_required_length << std::endl;
    return false;
  }

  const double length_tolerance = 0.20;
  double length_ratio = std::abs(diagonal1_length - diagonal2_length) / std::max(diagonal1_length, diagonal2_length);

  if (length_ratio > length_tolerance)
  {
    // std::cout << "  Cross REJECTED: diagonal length mismatch " << diagonal1_length << " vs " << diagonal2_length
    //           << " (ratio: " << length_ratio << ")" << std::endl;
    return false;
  }

  return true;
}

ToolPathPlanner::ConstPtr CrossSlicerRasterPlannerFactory::create() const
{
  auto planner = std::make_unique<CrossSlicerRasterPlanner>(direction_gen(), origin_gen());
  planner->setCrossSpacing(cross_spacing);
  planner->setLineSpacing(line_spacing);
  planner->setPointSpacing(point_spacing);
  planner->setMinHoleSize(min_hole_size);
  planner->setSearchRadius(search_radius);
  planner->setMinSegmentSize(min_segment_size);
  planner->generateRastersBidirectionally(bidirectional);

  return std::move(planner);
}

}  // namespace noether
