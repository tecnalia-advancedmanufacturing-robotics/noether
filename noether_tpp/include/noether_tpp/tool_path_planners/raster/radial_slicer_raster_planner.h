/**
 * @file radial_slicer_raster_planner.h
 * @copyright Copyright (c) 2021, Southwest Research Institute
 *
 * @par License
 * Software License Agreement (Apache License)
 * @par
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * @par
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once

#include <noether_tpp/tool_path_planners/raster/raster_planner.h>
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

namespace
{
struct RadialRasterConstructData
{
  std::vector<vtkSmartPointer<vtkPolyData>> raster_segments;
  std::vector<double> segment_lengths;
};

}  // namespace

namespace noether
{
/**
 * @ingroup radial_raster_planners
 * @brief An implementation of the Raster Planner using a series of cutting planes to generate radial cuts.
 * @details This implementation works best on approximately planar parts.
 * The direction generator defines the direction of the raster cut.
 * The cut normal (i.e., the raster step direction) is defined by the cross product of the cut direction and the
 * smallest principal axis of the mesh.
 */

class RadialSlicerRasterPlanner : public RasterPlanner
{
public:
  RadialSlicerRasterPlanner(DirectionGenerator::ConstPtr dir_gen, OriginGenerator::ConstPtr origin_gen);

  void setSearchRadius(const double search_radius);
  void setMinSegmentSize(const double min_segment_size);
  void generateRastersBidirectionally(const bool bidirectional);
  void setRadialNumCuts(const int num_radial_cuts);
  void setCutSymmetry(const bool symmetric_cuts);

protected:
  /**
   * @brief Implementation of the tool path planning capability
   * @param mesh
   * @return
   */
  ToolPaths planImpl(const pcl::PolygonMesh& mesh) const;

  /** @brief Flag indicating whether rasters should be generated in the direction of both the cut normal and its
   * negation */
  bool bidirectional_ = true;
  /** @brief Minimum length of valid segment (m) */
  double min_segment_size_;
  /** @brief Search radius for calculating normals (m) */
  double search_radius_;
  /** @brief Number of radial cuts */
  int num_radial_cuts_;
  /** @brief Symmetry of cuts */
  bool symmetric_cuts_;

  vtkSmartPointer<vtkPolyData> updateMesh(const pcl::PolygonMesh& mesh) const;
  void computeCuttingPlaneParameters(const pcl::PolygonMesh& mesh, Eigen::Vector3d& mesh_normal,
                                     Eigen::Matrix3d& pca_vecs, Eigen::Vector3d& centroid,
                                     Eigen::Vector3d& cut_direction, Eigen::Vector3d& cut_normal,
                                     Eigen::Vector3d& cut_origin) const;
  std::pair<double, double> calculateCuttingRange(const Eigen::Matrix3d& pca_vecs, const Eigen::Vector3d& centroid,
                                                  const Eigen::Vector3d& cut_origin,
                                                  const Eigen::Vector3d& cut_normal) const;
  RadialRasterConstructData processRasterSlice(const Eigen::Vector3d& cut_origin,
                                               vtkSmartPointer<vtkPolyData> raster_lines,
                                               vtkSmartPointer<vtkPolyData> mesh_data,
                                               const Eigen::Vector3d& cut_direction,
                                               vtkSmartPointer<vtkKdTreePointLocator> kd_tree,
                                               vtkSmartPointer<vtkCellLocator> cell_locator,
                                               std::size_t slice_index) const;
  std::vector<std::vector<vtkIdType>> extractRasterSegmentIds(vtkSmartPointer<vtkPolyData> raster_lines) const;
  vtkSmartPointer<vtkPoints> processAndAlignPoints(const std::vector<vtkIdType>& point_ids,
                                                   vtkSmartPointer<vtkPolyData> raster_lines,
                                                   const Eigen::Vector3d& cut_direction) const;
};

struct RadialSlicerRasterPlannerFactory : public RasterPlannerFactory
{
  bool bidirectional;
  double min_segment_size;
  double search_radius;
  int num_radial_cuts;
  bool symmetric_cuts;

  ToolPathPlanner::ConstPtr create() const override;
};

}  // namespace noether
