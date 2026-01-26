#include <noether_tpp/tool_path_modifiers/grid_raster_organization_modifier.h>
#include <noether_tpp/utils.h>

#include <numeric>

namespace noether
{
ToolPaths GridRasterOrganizationModifier::modify(ToolPaths tool_paths) const
{
  if (tool_paths.empty())
  {
    return tool_paths;
  }

  // FIXME: Sort the waypoints within each tool path segment by their x-axis value

  return tool_paths;
}

}  // namespace noether
