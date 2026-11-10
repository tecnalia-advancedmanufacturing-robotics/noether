#include <noether_gui/widgets/tool_path_planners/raster/direction_generators/fixed_direction_generator_widget.h>
#include "../../../ui_vector3d_editor_widget.h"
#include <noether_gui/utils.h>

#include <noether_tpp/tool_path_planners/raster/direction_generators/fixed_direction_generator.h>
#include <yaml-cpp/yaml.h>

#include <QVBoxLayout>

namespace
{
const std::string CUT_DIRECTION_KEY = "cut_direction";
const std::string CUT_NORMAL_KEY = "cut_normal";

void setVector(Ui::Vector3dEditor* ui, const Eigen::Vector3d& vector)
{
  ui->double_spin_box_x->setValue(vector.x());
  ui->double_spin_box_y->setValue(vector.y());
  ui->double_spin_box_z->setValue(vector.z());
}

Eigen::Vector3d readVector(const Ui::Vector3dEditor* ui)
{
  Eigen::Vector3d vector(
      ui->double_spin_box_x->value(), ui->double_spin_box_y->value(), ui->double_spin_box_z->value());
  vector.normalize();
  return vector;
}

Eigen::Vector3d loadVector(const YAML::Node& node)
{
  Eigen::Vector3d vector(noether::getEntry<double>(node, "x"),
                         noether::getEntry<double>(node, "y"),
                         noether::getEntry<double>(node, "z"));
  vector.normalize();
  return vector;
}

void saveVector(YAML::Node& node, const Ui::Vector3dEditor* ui)
{
  node["x"] = ui->double_spin_box_x->value();
  node["y"] = ui->double_spin_box_y->value();
  node["z"] = ui->double_spin_box_z->value();
}
}  // namespace

namespace noether
{
FixedDirectionGeneratorWidget::FixedDirectionGeneratorWidget(QWidget* parent)
  : DirectionGeneratorWidget(parent)
  , cut_direction_ui_(new Ui::Vector3dEditor())
  , cut_normal_ui_(new Ui::Vector3dEditor())
{
  auto* layout = new QVBoxLayout(this);

  auto* direction_page = new QWidget(this);
  cut_direction_ui_->setupUi(direction_page);
  cut_direction_ui_->group_box->setTitle("Cut Direction");
  cut_direction_ui_->double_spin_box_x->setValue(1.0);
  layout->addWidget(direction_page);

  auto* normal_page = new QWidget(this);
  cut_normal_ui_->setupUi(normal_page);
  cut_normal_ui_->group_box->setTitle("Mesh Normal");
  cut_normal_ui_->double_spin_box_z->setValue(1.0);
  layout->addWidget(normal_page);
}

void FixedDirectionGeneratorWidget::configure(const YAML::Node& config)
{
  setVector(cut_direction_ui_, loadVector(config[CUT_DIRECTION_KEY]));
  setVector(cut_normal_ui_, loadVector(config[CUT_NORMAL_KEY]));
}

void FixedDirectionGeneratorWidget::save(YAML::Node& config) const
{
  YAML::Node cut_direction;
  saveVector(cut_direction, cut_direction_ui_);
  config[CUT_DIRECTION_KEY] = cut_direction;

  YAML::Node cut_normal;
  saveVector(cut_normal, cut_normal_ui_);
  config[CUT_NORMAL_KEY] = cut_normal;
}

DirectionGenerator::ConstPtr FixedDirectionGeneratorWidget::create() const
{
  return std::make_unique<FixedDirectionGenerator>(readVector(cut_direction_ui_), readVector(cut_normal_ui_));
}

}  // namespace noether
