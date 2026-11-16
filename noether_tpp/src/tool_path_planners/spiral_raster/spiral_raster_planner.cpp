#include <noether_tpp/tool_path_planners/spiral_raster/spiral_raster_planner.h>
#include <noether_tpp/tool_path_modifiers/direction_of_travel_orientation_modifier.h>

#include <utility>

namespace noether
{
SpiralRasterPlanner::SpiralRasterPlanner(SpiralDirectionGenerator::ConstPtr dir_gen) : dir_gen_(std::move(dir_gen))
{
}

ToolPaths SpiralRasterPlanner::plan(const pcl::PolygonMesh& mesh) const
{
  ToolPaths tool_paths = planImpl(mesh);
  if (tool_paths.empty())
    return tool_paths;

  // Keep the z-axis from the planner (surface normal) and align x with travel along the spiral.
  DirectionOfTravelOrientationModifier dot_orientation;
  return dot_orientation.modify(tool_paths);
}

void SpiralRasterPlanner::setPointSpacing(const double point_spacing) { point_spacing_ = point_spacing; }
void SpiralRasterPlanner::setLineSpacing(const double line_spacing) { line_spacing_ = line_spacing; }
void SpiralRasterPlanner::setMinHoleSize(const double min_hole_size) { min_hole_size_ = min_hole_size; }

}  // namespace noether
