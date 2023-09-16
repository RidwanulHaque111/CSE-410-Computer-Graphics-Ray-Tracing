#include "1705111_vector.hpp"

class Ray
{
public:
    VectorPoint start;
    VectorPoint dir;
    Ray(VectorPoint, VectorPoint);
    VectorPoint getStartPoint();
    VectorPoint getDirPoint();
};

Ray::Ray(VectorPoint from, VectorPoint to)
{
    start = from;
    dir = to;
    dir.normalize();
}

VectorPoint Ray::getStartPoint()
{
    return start;
}

VectorPoint Ray::getDirPoint()
{
    return dir;
}