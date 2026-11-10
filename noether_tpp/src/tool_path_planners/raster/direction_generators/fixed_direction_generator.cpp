#include <noether_tpp/tool_path_planners/raster/direction_generators/fixed_direction_generator.h>

namespace noether
{
FixedDirectionGenerator::FixedDirectionGenerator(const Eigen::Vector3d& cut_direction, const Eigen::Vector3d& cut_normal) : cut_direction_(cut_direction.normalized()), cut_normal_(cut_normal.normalized())
{
}

Eigen::Vector3d FixedDirectionGenerator::generate_cut_direction(const pcl::PolygonMesh&) const { return cut_direction_; }

Eigen::Vector3d FixedDirectionGenerator::generate_mesh_normal(const pcl::PolygonMesh&) const { return cut_normal_; }

}  // namespace noether
