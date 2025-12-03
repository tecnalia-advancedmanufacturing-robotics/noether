#include <noether_tpp/tool_path_modifiers/cross_raster_organization_modifier.h>
#include <noether_tpp/utils.h>

#include <numeric>

namespace noether
{
ToolPaths CrossRasterOrganizationModifier::modify(ToolPaths tool_paths) const
{
  if (tool_paths.empty())
  {
    return tool_paths;
  }

  // Sort the waypoints within each diagonal segment by their distance along the diagonal direction
  for (ToolPath& cross : tool_paths)
  {
    for (ToolPathSegment& diagonal : cross)
    {
      if (diagonal.size() < 2)
        continue;

      const Eigen::Vector3d diagonal_start = diagonal.front().translation();
      const Eigen::Vector3d diagonal_end = diagonal.back().translation();
      const Eigen::Vector3d diagonal_direction = (diagonal_end - diagonal_start).normalized();

      std::sort(diagonal.begin(), diagonal.end(),
                [&diagonal_start, &diagonal_direction](const ToolPathWaypoint& a, const ToolPathWaypoint& b) {
                  Eigen::Vector3d diff_a = a.translation() - diagonal_start;
                  Eigen::Vector3d diff_b = b.translation() - diagonal_start;
                  return diff_a.dot(diagonal_direction) < diff_b.dot(diagonal_direction);
                });
    }
  }

  // Sort diagonals within each cross by the angle their direction makes in the XY plane
  for (ToolPath& cross : tool_paths)
  {
    if (cross.size() < 2)
      continue;

    std::sort(cross.begin(), cross.end(), [](const ToolPathSegment& a, const ToolPathSegment& b) {
      if (a.size() < 2 || b.size() < 2)
        return a.size() >= 2;

      Eigen::Vector3d dir_a = (a.back().translation() - a.front().translation()).normalized();
      Eigen::Vector3d dir_b = (b.back().translation() - b.front().translation()).normalized();
      double angle_a = std::atan2(dir_a.y(), dir_a.x());
      double angle_b = std::atan2(dir_b.y(), dir_b.x());

      return angle_a < angle_b;
    });
  }

  // Sort crosses by the position of a reference point
  if (tool_paths.size() > 1)
  {
    Eigen::Vector3d reference_point = tool_paths.at(0).at(0).front().translation();
    Eigen::Vector3d raster_direction = Eigen::Vector3d::UnitX();  // FIXME: correct default?
    if (tool_paths.size() >= 2 && !tool_paths.at(1).empty() && !tool_paths.at(1).at(0).empty())
    {
      Eigen::Vector3d second_cross_point = tool_paths.at(1).at(0).front().translation();
      raster_direction = (second_cross_point - reference_point).normalized();
    }
    else if (tool_paths.size() >= 3)
    {
      std::vector<Eigen::Vector3d> first_points;
      for (const auto& cross : tool_paths)
      {
        if (!cross.empty() && !cross.at(0).empty())
        {
          first_points.push_back(cross.at(0).front().translation());
        }
      }

      if (first_points.size() >= 2)
      {
        Eigen::Vector3d avg_direction = Eigen::Vector3d::Zero();
        for (size_t i = 1; i < first_points.size(); ++i)
        {
          avg_direction += (first_points[i] - first_points[i - 1]).normalized();
        }
        raster_direction = avg_direction.normalized();
      }
    }

    // Sort crosses by projecting the first point of diagonal1 onto the raster direction
    std::sort(tool_paths.begin(), tool_paths.end(),
              [&reference_point, &raster_direction](const ToolPath& cross_a, const ToolPath& cross_b) {
                if (cross_a.empty() || cross_a.at(0).empty())
                  return false;
                if (cross_b.empty() || cross_b.at(0).empty())
                  return true;

                Eigen::Vector3d point_a = cross_a.at(0).front().translation();
                Eigen::Vector3d point_b = cross_b.at(0).front().translation();
                double proj_a = (point_a - reference_point).dot(raster_direction);
                double proj_b = (point_b - reference_point).dot(raster_direction);
                return proj_a < proj_b;
              });
  }

  return tool_paths;
}

}  // namespace noether
