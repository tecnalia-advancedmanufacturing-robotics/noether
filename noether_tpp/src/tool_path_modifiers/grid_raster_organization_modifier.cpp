#include <noether_tpp/tool_path_modifiers/grid_raster_organization_modifier.h>
#include <noether_tpp/utils.h>

#include <numeric>

namespace noether
{
ToolPaths GridRasterOrganizationModifier::modify(ToolPaths tool_paths) const
{
  // Find the index correspondig to the change of direction in the grid
  std::vector<Eigen::Vector3d> path_directions;
  for (size_t i = 0; i < tool_paths.size(); ++i)
  {
    path_directions.push_back(estimateToolPathDirection(tool_paths.at(i)));
  }

  int index_change_orientation = -1;
  Eigen::Vector3d reference_segment_dir_1 = path_directions[0];
  Eigen::Vector3d reference_segment_dir_2;

  for (size_t i = 1; i < path_directions.size(); ++i)
  {
    // Check if direction has changed
    if (std::abs(path_directions[i].dot(reference_segment_dir_1)) < 0.5)
    {
      index_change_orientation = i;
      reference_segment_dir_2 = path_directions[i];
      break;
    }
  }

  int i = 0;
  for (ToolPath& tool_path : tool_paths)
  {
    const Eigen::Vector3d reference_segment_dir = path_directions[i];

    // Sort the waypoints within each tool path segment by their distance along the reference direction
    for (ToolPathSegment& segment : tool_path)
    {
      std::sort(segment.begin(), segment.end(),
                [&segment, &reference_segment_dir](const ToolPathWaypoint& a, const ToolPathWaypoint& b) {
                  Eigen::Vector3d diff_from_start_b = b.translation() - segment.at(0).translation();
                  Eigen::Vector3d diff_from_start_a = a.translation() - segment.at(0).translation();
                  return diff_from_start_a.dot(reference_segment_dir) < diff_from_start_b.dot(reference_segment_dir);
                });
    }

    // Sort the tool path segments within each tool path by the distance of their first waypoints along the reference
    // direction
    std::sort(tool_path.begin(), tool_path.end(),
              [&tool_path, &reference_segment_dir](const ToolPathSegment& a, const ToolPathSegment& b) {
                Eigen::Vector3d diff_from_start_b = b.at(0).translation() - tool_path.at(0).at(0).translation();
                Eigen::Vector3d diff_from_start_a = a.at(0).translation() - tool_path.at(0).at(0).translation();
                return diff_from_start_a.dot(reference_segment_dir) < diff_from_start_b.dot(reference_segment_dir);
              });

    i++;
  }

  if (index_change_orientation == -1)
  {
    // Single orientation - Same as raster_organization_modifier
    // Sort the tool paths by their distance along a vector that is perpendicular to the reference direction of travel
    const Eigen::Vector3d reference_tool_paths_dir = estimateRasterDirection(tool_paths, reference_segment_dir_1);
    const Eigen::Isometry3d first_wp = tool_paths.at(0).at(0).at(0);

    std::sort(tool_paths.begin(), tool_paths.end(),
              [&reference_tool_paths_dir, &first_wp](const ToolPath& a, const ToolPath& b) {
                Eigen::Vector3d diff_from_start_a = a.at(0).at(0).translation() - first_wp.translation();
                Eigen::Vector3d diff_from_start_b = b.at(0).at(0).translation() - first_wp.translation();
                return diff_from_start_a.dot(reference_tool_paths_dir) <
                       diff_from_start_b.dot(reference_tool_paths_dir);
              });
  }
  else
  {
    // Two orientations - Grid Raster
    // Sort the tool paths by their distance along a vector that is perpendicular to the reference direction of travel
    // Sort each group separately and then combine
    ToolPaths group1(tool_paths.begin(), tool_paths.begin() + index_change_orientation);
    ToolPaths group2(tool_paths.begin() + index_change_orientation, tool_paths.end());

    // Sort first group
    const Eigen::Vector3d reference_tool_paths_dir_1 = estimateRasterDirection(group1, reference_segment_dir_1);
    const Eigen::Isometry3d first_wp_1 = group1.at(0).at(0).at(0);

    std::sort(
        group1.begin(), group1.end(), [&reference_tool_paths_dir_1, &first_wp_1](const ToolPath& a, const ToolPath& b) {
          Eigen::Vector3d diff_from_start_a = a.at(0).at(0).translation() - first_wp_1.translation();
          Eigen::Vector3d diff_from_start_b = b.at(0).at(0).translation() - first_wp_1.translation();
          return diff_from_start_a.dot(reference_tool_paths_dir_1) < diff_from_start_b.dot(reference_tool_paths_dir_1);
        });

    // Sort second group
    const Eigen::Vector3d reference_tool_paths_dir_2 = estimateRasterDirection(group2, reference_segment_dir_2);
    const Eigen::Isometry3d first_wp_2 = group2.at(0).at(0).at(0);

    std::sort(
        group2.begin(), group2.end(), [&reference_tool_paths_dir_2, &first_wp_2](const ToolPath& a, const ToolPath& b) {
          Eigen::Vector3d diff_from_start_a = a.at(0).at(0).translation() - first_wp_2.translation();
          Eigen::Vector3d diff_from_start_b = b.at(0).at(0).translation() - first_wp_2.translation();
          return diff_from_start_a.dot(reference_tool_paths_dir_2) < diff_from_start_b.dot(reference_tool_paths_dir_2);
        });

    // Combine
    tool_paths.clear();
    tool_paths.insert(tool_paths.end(), group1.begin(), group1.end());
    tool_paths.insert(tool_paths.end(), group2.begin(), group2.end());
  }

  return tool_paths;
}

}  // namespace noether
