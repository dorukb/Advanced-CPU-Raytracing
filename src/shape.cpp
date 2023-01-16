#include "shape.hpp"

using namespace DorkTracer;

int DorkTracer::Shape::GetMaterial(){
    return this->material_id;
}
void DorkTracer::Shape::SetMaterial(int matId){
    this->material_id = matId;
}