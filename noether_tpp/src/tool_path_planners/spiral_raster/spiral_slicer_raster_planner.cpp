#include <noether_tpp/tool_path_planners/spiral_raster/spiral_slicer_raster_planner.h>
#include <noether_tpp/utils.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include <pcl/common/common.h>
#include <pcl/common/pca.h>
#include <pcl/conversions.h>
#include <pcl/surface/vtk_smoothing/vtk_utils.h>

#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>

namespace
{
constexpr double kEpsilon = 1e-6;
// Near the center the radius is smaller than the point spacing; keep at least 8 samples per turn there.
constexpr double kMaxThetaStep = M_PI / 4.0;

Eigen::Matrix3d frame(const Eigen::Vector3d& x_axis, const Eigen::Vector3d& y_axis, const Eigen::Vector3d& z_axis)
{
  Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
  rotation.col(0) = x_axis.normalized();
  rotation.col(1) = y_axis.normalized();
  rotation.col(2) = z_axis.normalized();
  return rotation;
}

void principalAxes(const pcl::PolygonMesh& mesh, Eigen::Vector3d& centroid, Eigen::Matrix3d& axes)
{
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
  pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

  pcl::PCA<pcl::PointXYZ> pca;
  pca.setInputCloud(cloud);

  pcl::PointCloud<pcl::PointXYZ> projected;
  pca.project(*cloud, projected);
  pcl::PointXYZ min_point;
  pcl::PointXYZ max_point;
  pcl::getMinMax3D(projected, min_point, max_point);
  const Eigen::Array3f scales = max_point.getArray3fMap() - min_point.getArray3fMap();

  centroid = pca.getMean().head<3>().cast<double>();
  axes = (pca.getEigenVectors().array().rowwise() * scales.transpose()).cast<double>();
  if (axes.col(0).sum() < 0.0)
    axes.col(0) *= -1.0;
  if (axes.col(1).sum() < 0.0)
    axes.col(1) *= -1.0;
  axes.col(2) = axes.col(0).cross(axes.col(1)).normalized() * scales(2);
}

double spiralReach(const Eigen::Matrix3d& axes, const Eigen::Vector3d& centroid, const Eigen::Vector3d& center,
                   const Eigen::Vector3d& x_axis, const Eigen::Vector3d& y_axis)
{
  Eigen::MatrixXi corner_mask(3, 8);
  // clang-format off
  corner_mask <<
      1, -1, 1, -1, 1, -1, 1, -1,
      1, 1, -1, -1, 1, 1, -1, -1,
      1, 1, 1, 1, -1, -1, -1, -1;
  // clang-format on

  Eigen::MatrixXd corners = (axes / 2.0) * corner_mask.cast<double>();
  corners.colwise() += centroid;

  double reach = 0.0;
  for (Eigen::Index col = 0; col < corners.cols(); ++col)
  {
    const Eigen::Vector3d offset = corners.col(col) - center;
    reach = std::max(reach, std::hypot(offset.dot(x_axis), offset.dot(y_axis)));
  }
  return reach;
}

bool projectOntoMesh(vtkCellLocator* locator, const Eigen::Vector3d& plane_point, const Eigen::Vector3d& mesh_normal,
                     double ray_length, Eigen::Vector3d& hit, vtkIdType& cell_id)
{
  const double tolerance = 1e-6;
  const double approach = 1e-5;
  double best = std::numeric_limits<double>::max();
  bool found = false;

  auto consider = [&](const Eigen::Vector3d& start, const Eigen::Vector3d& end) {
    double t = 0.0;
    double intersection[3];
    double pcoords[3];
    int sub_id = 0;
    vtkIdType candidate_cell = -1;
    if (locator->IntersectWithLine(start.data(), end.data(), tolerance, t, intersection, pcoords, sub_id,
                                   candidate_cell) == 0)
    {
      return;
    }

    const Eigen::Vector3d candidate(intersection[0], intersection[1], intersection[2]);
    const double distance = (candidate - plane_point).norm();
    if (distance < best)
    {
      best = distance;
      hit = candidate;
      cell_id = candidate_cell;
      found = true;
    }
  };

  const Eigen::Vector3d normal = mesh_normal.normalized();
  consider(plane_point - normal * approach, plane_point + normal * ray_length);
  consider(plane_point + normal * approach, plane_point - normal * ray_length);
  return found;
}

Eigen::Vector3d surfaceNormal(double search_radius, vtkPolyData* mesh, vtkKdTreePointLocator* kd_tree,
                              vtkCellLocator* cell_locator, const Eigen::Vector3d& point, vtkIdType cell_id,
                              const Eigen::Vector3d& reference_normal)
{
  Eigen::Vector3d normal = Eigen::Vector3d::Zero();
  vtkDataArray* vertex_normals = mesh->GetPointData()->GetNormals();

  if (search_radius > 0.0 && vertex_normals != nullptr)
  {
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    kd_tree->FindPointsWithinRadius(search_radius, point.data(), ids);
    if (ids->GetNumberOfIds() < 1)
      kd_tree->FindClosestNPoints(1, point.data(), ids);

    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i)
    {
      const vtkIdType id = ids->GetId(i);
      if (id < 0)
        continue;
      Eigen::Vector3d sample;
      vertex_normals->GetTuple(id, sample.data());
      normal += sample.normalized();
    }
  }

  if (normal.norm() <= kEpsilon && cell_id >= 0 && mesh->GetCellData()->GetNormals() != nullptr)
    mesh->GetCellData()->GetNormals()->GetTuple(cell_id, normal.data());

  if (normal.norm() <= kEpsilon)
  {
    Eigen::Vector3d closest;
    vtkIdType closest_cell = -1;
    int sub_id = 0;
    double distance2 = 0.0;
    cell_locator->FindClosestPoint(point.data(), closest.data(), closest_cell, sub_id, distance2);
    if (closest_cell >= 0 && mesh->GetCellData()->GetNormals() != nullptr)
      mesh->GetCellData()->GetNormals()->GetTuple(closest_cell, normal.data());
  }

  if (normal.norm() <= kEpsilon)
    throw std::runtime_error("Could not estimate a surface normal for a spiral waypoint");

  normal.normalize();
  if (normal.dot(reference_normal) < 0.0)
    normal *= -1.0;
  return normal;
}

Eigen::Vector3d averageVertexNormal(vtkPolyData* mesh)
{
  Eigen::Vector3d average = Eigen::Vector3d::Zero();
  vtkDataArray* vertex_normals = mesh->GetPointData()->GetNormals();
  for (vtkIdType id = 0; id < vertex_normals->GetNumberOfTuples(); ++id)
  {
    Eigen::Vector3d sample;
    vertex_normals->GetTuple(id, sample.data());
    average += sample;
  }
  return average;
}

vtkSmartPointer<vtkPolyData> meshWithNormals(const pcl::PolygonMesh& mesh)
{
  if (!noether::hasNormals(mesh))
  {
    throw std::runtime_error(
        "The input mesh does not have vertex normals. Use a MeshModifier to generate them before planning a spiral.");
  }

  vtkSmartPointer<vtkPolyData> mesh_data = vtkSmartPointer<vtkPolyData>::New();
  pcl::VTKUtils::mesh2vtk(mesh, mesh_data);
  mesh_data->BuildLinks();
  mesh_data->BuildCells();

  if (mesh_data->GetPointData()->GetNormals() == nullptr || mesh_data->GetCellData()->GetNormals() == nullptr)
  {
    vtkSmartPointer<vtkPolyDataNormals> normal_generator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_generator->SetInputData(mesh_data);
    normal_generator->ComputePointNormalsOn();
    normal_generator->SetComputeCellNormals(mesh_data->GetCellData()->GetNormals() == nullptr);
    normal_generator->SetFeatureAngle(M_PI_2);
    normal_generator->SetSplitting(true);
    normal_generator->SetConsistency(true);
    normal_generator->SetAutoOrientNormals(false);
    normal_generator->SetFlipNormals(false);
    normal_generator->SetNonManifoldTraversal(false);
    normal_generator->Update();

    if (mesh_data->GetPointData()->GetNormals() == nullptr)
      mesh_data->GetPointData()->SetNormals(normal_generator->GetOutput()->GetPointData()->GetNormals());
    if (mesh_data->GetCellData()->GetNormals() == nullptr)
      mesh_data->GetCellData()->SetNormals(normal_generator->GetOutput()->GetCellData()->GetNormals());
  }

  return mesh_data;
}

}  // namespace

namespace noether
{
SpiralSlicerRasterPlanner::SpiralSlicerRasterPlanner(SpiralDirectionGenerator::ConstPtr dir_gen)
  : SpiralRasterPlanner(std::move(dir_gen)), search_radius_(0.0), min_segment_size_(0.0), delta_x_(0.0), delta_y_(0.0)
{
}

void SpiralSlicerRasterPlanner::setSearchRadius(const double search_radius) { search_radius_ = search_radius; }

void SpiralSlicerRasterPlanner::setMinSegmentSize(const double min_segment_size)
{
  min_segment_size_ = min_segment_size;
}

void SpiralSlicerRasterPlanner::setCenterOffset(const double delta_x, const double delta_y)
{
  delta_x_ = delta_x;
  delta_y_ = delta_y;
}

ToolPaths SpiralSlicerRasterPlanner::planImpl(const pcl::PolygonMesh& mesh) const
{
  if (line_spacing_ <= kEpsilon || point_spacing_ <= kEpsilon)
    return {};

  vtkSmartPointer<vtkPolyData> mesh_data = meshWithNormals(mesh);

  vtkSmartPointer<vtkKdTreePointLocator> kd_tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
  kd_tree->SetDataSet(mesh_data);
  kd_tree->BuildLocator();

  vtkSmartPointer<vtkCellLocator> cell_locator = vtkSmartPointer<vtkCellLocator>::New();
  cell_locator->SetDataSet(mesh_data);
  cell_locator->BuildLocator();

  Eigen::Vector3d centroid;
  Eigen::Matrix3d axes;
  principalAxes(mesh, centroid, axes);

  // Generators such as the principal axis one return a normal with an arbitrary sign; orient it like the mesh.
  Eigen::Vector3d mesh_normal = dir_gen_->generate_mesh_normal(mesh).normalized();
  if (mesh_normal.dot(averageVertexNormal(mesh_data)) < 0.0)
    mesh_normal *= -1.0;
  Eigen::Vector3d x_axis = dir_gen_->generate_cut_direction(mesh);
  x_axis = (x_axis - mesh_normal * x_axis.dot(mesh_normal)).normalized();
  const Eigen::Vector3d y_axis = mesh_normal.cross(x_axis).normalized();

  const Eigen::Vector3d center = centroid + delta_x_ * x_axis + delta_y_ * y_axis;
  const double radius_limit = spiralReach(axes, centroid, center, x_axis, y_axis);
  const double pitch = line_spacing_ / (2.0 * M_PI);
  const double theta_limit = radius_limit / pitch;
  const double ray_length = std::max(axes.colwise().norm().maxCoeff(), radius_limit);
  const double gap_distance = std::max(min_hole_size_, 2.0 * point_spacing_);

  ToolPath spiral;
  ToolPathSegment segment;
  bool has_previous = false;
  Eigen::Vector3d previous_point = Eigen::Vector3d::Zero();
  double segment_length = 0.0;

  auto flush = [&]() {
    if (segment.size() > 1 && segment_length > min_segment_size_)
      spiral.push_back(segment);
    segment.clear();
    segment_length = 0.0;
    has_previous = false;
  };

  for (double theta = 0.0; theta <= theta_limit + kEpsilon;)
  {
    const double radius = pitch * theta;
    const Eigen::Vector3d plane_point =
        center + (radius * std::cos(theta)) * x_axis + (radius * std::sin(theta)) * y_axis;

    Eigen::Vector3d hit;
    vtkIdType cell_id = -1;
    if (projectOntoMesh(cell_locator, plane_point, mesh_normal, ray_length, hit, cell_id))
    {
      if (has_previous && (hit - previous_point).norm() > gap_distance)
        flush();

      const Eigen::Vector3d normal =
          surfaceNormal(search_radius_, mesh_data, kd_tree, cell_locator, hit, cell_id, mesh_normal);
      Eigen::Vector3d travel = x_axis;
      if (has_previous && (hit - previous_point).norm() > kEpsilon)
        travel = (hit - previous_point).normalized();
      Eigen::Vector3d y = normal.cross(travel);
      if (y.norm() <= kEpsilon)
        y = normal.cross(x_axis);
      y.normalize();
      travel = y.cross(normal).normalized();

      if (!segment.empty())
        segment_length += (hit - segment.back().translation()).norm();

      segment.push_back(Eigen::Translation3d(hit) * Eigen::AngleAxisd(frame(travel, y, normal)));
      previous_point = hit;
      has_previous = true;
    }
    else if (has_previous)
    {
      flush();
    }

    const double arc_scale = std::sqrt(radius * radius + pitch * pitch);
    theta += std::min(point_spacing_ / arc_scale, kMaxThetaStep);
  }
  flush();

  if (spiral.empty())
    return {};

  ToolPaths tool_paths;
  tool_paths.push_back(std::move(spiral));
  return tool_paths;
}

ToolPathPlanner::ConstPtr SpiralSlicerRasterPlannerFactory::create() const
{
  auto planner = std::make_unique<SpiralSlicerRasterPlanner>(direction_gen());
  planner->setLineSpacing(line_spacing);
  planner->setPointSpacing(point_spacing);
  planner->setMinHoleSize(min_hole_size);
  planner->setSearchRadius(search_radius);
  planner->setMinSegmentSize(min_segment_size);
  planner->setCenterOffset(delta_x, delta_y);
  return planner;
}

}  // namespace noether
