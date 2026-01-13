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

  return tool_paths;
}

}  // namespace noether
