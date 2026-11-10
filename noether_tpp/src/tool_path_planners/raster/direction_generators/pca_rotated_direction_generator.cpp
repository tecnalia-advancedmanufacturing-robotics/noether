#include <noether_tpp/tool_path_planners/raster/direction_generators/pca_rotated_direction_generator.h>

#include <pcl/common/pca.h>

namespace noether
{
PCARotatedDirectionGenerator::PCARotatedDirectionGenerator(DirectionGenerator::ConstPtr&& dir_gen,
                                                           double rotation_angle)
  : dir_gen_(std::move(dir_gen)), rotation_angle_(rotation_angle)
{
}

Eigen::Vector3d PCARotatedDirectionGenerator::generate_cut_direction(const pcl::PolygonMesh& mesh) const
{
  Eigen::Vector3d nominal_dir = dir_gen_->generate_cut_direction(mesh);

  auto vertices = pcl::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
  pcl::fromPCLPointCloud2(mesh.cloud, *vertices);
  pcl::PCA<pcl::PointXYZ> pca;
  pca.setInputCloud(vertices);
  Eigen::Matrix3d pca_vecs = pca.getEigenVectors().cast<double>();

  // Rotate the nominal direction about the smallest principal axis by the specified angle
  return Eigen::AngleAxisd(rotation_angle_, pca_vecs.col(2).normalized()) * nominal_dir;
}

Eigen::Vector3d PCARotatedDirectionGenerator::generate_mesh_normal(const pcl::PolygonMesh& mesh) const
{
  return dir_gen_->generate_mesh_normal(mesh);
}

}  // namespace noether
