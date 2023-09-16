#include <bits/stdc++.h>
using namespace std;

#define pi 3.14159265358979323846

#define ZERO 1e-8
#define INF 1e12


double rad(double deg)
{
    double converted = (deg * pi) / 180.0;
    return converted;
}

// vector class ==============================================
class VectorPoint
{
public:
    double x, y, z;

    VectorPoint() : x(0.0), y(0.0), z(0.0) {}

    VectorPoint(double x, double y, double z) : x(x), y(y), z(z) {}

    friend ostream &operator<<(ostream &out, const VectorPoint &vect)
    {
        out << "point is: (";
        out << vect.x << "," << vect.y << "," << vect.z;
        out << ")" << endl;
        return out;
    }

    friend istream &operator>>(istream &ins, VectorPoint &vect)
    {
        ins >> vect.x >> vect.y >> vect.z;
        return ins;
    }

    bool operator==(const VectorPoint& other) const
    {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }


    //add vectors
    VectorPoint operator+(const VectorPoint &vect) const
    {
        double x_ = this->x + vect.x;
        double y_ = this->y + vect.y;
        double z_ = this->z + vect.z;
        return VectorPoint(x_, y_, z_);
    }

    //add constant
    VectorPoint operator+(double val) const
    {
        double x_ = this->x + val;
        double y_ = this->y + val;
        double z_ = this->z + val;
        return VectorPoint(x_, y_, z_);
    }

    //VectorPoint * VectorPoint operator overloading
    VectorPoint operator*(const VectorPoint &vect) const
    {
        double x_ = this->x * vect.x;
        double y_ = this->y * vect.y;
        double z_ = this->z * vect.z;
        return VectorPoint(x_, y_, z_);
    }

    //sub vectors
    VectorPoint operator-(const VectorPoint &vect) const
    {
        double x_ = this->x - vect.x;
        double y_ = this->y - vect.y;
        double z_ = this->z - vect.z;
        return VectorPoint(x_, y_, z_);
    }

    //sub constant
    VectorPoint operator-(double val) const
    {
        double x_ = this->x - val;
        double y_ = this->y - val;
        double z_ = this->z - val;
        return VectorPoint(x_, y_, z_);
    }

    //negate vectors
    VectorPoint operator-() const
    {
        double x_ = -this->x;
        double y_ = -this->y;
        double z_ = -this->z;
        return VectorPoint(x_, y_, z_);
    }

    //mul constant
    VectorPoint operator*(double val) const
    {
        double x_ = this->x * val;
        double y_ = this->y * val;
        double z_ = this->z * val;
        return VectorPoint(x_, y_, z_);
    }

    //divide constant
    VectorPoint operator/(double val) const
    {
        double x_ = this->x / val;
        double y_ = this->y / val;
        double z_ = this->z / val;
        return VectorPoint(x_, y_, z_);
    }

    //dot product
    double operator^(const VectorPoint &vect) const
    {
        double x_ = this->x * vect.x;
        x_ += this->y * vect.y;
        x_ += this->z * vect.z;
        return x_;
    }

    //cross product
    VectorPoint cross(VectorPoint vect)
    {
        double x_  = this->y * vect.z - this->z * vect.y;
        double y_  = this->z * vect.x - this->x * vect.z;
        double z_  = this->x * vect.y - this->y * vect.x;
        return VectorPoint(x_, y_, z_);
    }

    //vector length
    double get_length()
    {
        double x_ = this->x * this->x;
        double y_ = this->y * this->y;
        double z_ = this->z * this->z;
        return sqrt(x_+y_+z_);
    }

    double magnitude()
    {
        return sqrt(x * x + y * y + z * z);
    }

    //normalize a vector
    void normalize()
    {
        double mag = magnitude();
        if (mag > 0)
        {
            x /= mag;
            y /= mag;
            z /= mag;
        }
    }

};


//utils from first offline =====================
VectorPoint get_VectorPoint(double x, double y, double z)
{
    VectorPoint p1;
    p1.x = x;
    p1.y = y;
    p1.z = z;
    return p1;
}

VectorPoint vect_add(VectorPoint p1, VectorPoint p2)
{
    p1.x = p1.x + p2.x;
    p1.y = p1.y + p2.y;
    p1.z = p1.z + p2.z;
    return p1;
}

VectorPoint vect_scale(VectorPoint p1, double val)
{
    p1.x = p1.x * val;
    p1.y = p1.y * val;
    p1.z = p1.z * val;
    return p1;
}

VectorPoint unit_vector(VectorPoint p1)
{
    double mod = 1/sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
    return vect_scale(p1, mod);
}

VectorPoint move_along_vect(VectorPoint p, VectorPoint vect, double units)
{
    VectorPoint dir = vect_add(vect, vect_scale(p, -1)); //direction of move
    dir = vect_scale(unit_vector(dir), units); //units of movement
    p = vect_add(p, dir);
    return p;
}

VectorPoint move_along_unit_vect(VectorPoint p, VectorPoint uvect, double units)
{
    VectorPoint dir = vect_scale(uvect, units); //units of movement
    p = vect_add(p, dir);
    return p;
}

VectorPoint cross_product(VectorPoint a, VectorPoint b)
{
    VectorPoint p;
    p.x = a.y*b.z - a.z*b.y;
    p.y = a.z*b.x - a.x*b.z;
    p.z = a.x*b.y - a.y*b.x;
    return p;
}



void axis_rotation(VectorPoint *l, VectorPoint *u, VectorPoint *r, double degrees)
{
    double rad = degrees * pi / 180.0;
    *l = vect_add(vect_scale(*l, cos(rad)), vect_scale(*u, sin(rad)));
    *u = cross_product(*r, *l);
}

VectorPoint rotate_about_VectorPoint(VectorPoint d, VectorPoint o, double degrees)
{
    double rad = degrees * pi / 180.0;
    d.x -= o.x;
    d.y -= o.y;

    double temp = d.x*cos(rad) - d.y*sin(rad);
    d.y = d.x*sin(rad) + d.y*cos(rad);
    d.x = temp;

    d.x += o.x;
    d.y += o.y;
    return d;
}

//quad equation solution
VectorPoint solve_quad(double a, double b, double c)
{
    double discrim = b*b-4*a*c;
    if(discrim < 0) return VectorPoint(-1, -1, -1);
    double sol1 = (-b-sqrt(discrim))/(2*a);
    double sol2 = (-b+sqrt(discrim))/(2*a);
    return VectorPoint(sol1, sol2, 0);
}


double calc_determinant(double mat[3][3])
{
    double ans;
    ans = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
          - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
          + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return ans;
}

// hard coded triple linear equations using cramer's rule
VectorPoint solve_tri_linears(double a, double b, double c, double d,
                           double e, double f, double g, double h,
                           double i, double j, double k, double l)
{
    double coeff[3][4] =
    {
        { a, b, c, d },
        { e, f, g, h },
        { i, j, k, l },
    };

    double det[3][3] =
    {
        { coeff[0][0], coeff[0][1], coeff[0][2] },
        { coeff[1][0], coeff[1][1], coeff[1][2] },
        { coeff[2][0], coeff[2][1], coeff[2][2] },
    };
    double d1[3][3] =
    {
        { coeff[0][3], coeff[0][1], coeff[0][2] },
        { coeff[1][3], coeff[1][1], coeff[1][2] },
        { coeff[2][3], coeff[2][1], coeff[2][2] },
    };
    double d2[3][3] =
    {
        { coeff[0][0], coeff[0][3], coeff[0][2] },
        { coeff[1][0], coeff[1][3], coeff[1][2] },
        { coeff[2][0], coeff[2][3], coeff[2][2] },
    };
    double d3[3][3] =
    {
        { coeff[0][0], coeff[0][1], coeff[0][3] },
        { coeff[1][0], coeff[1][1], coeff[1][3] },
        { coeff[2][0], coeff[2][1], coeff[2][3] },
    };

    double D = calc_determinant(det);
    double D1 = calc_determinant(d1);
    double D2 = calc_determinant(d2);
    double D3 = calc_determinant(d3);

    // if sol exists
    if (D != 0)
    {
        return VectorPoint(D1/D, D2/D, D3/D);
    }

    //return invalid k1, k2, tmin
    return VectorPoint(10, 10, -1.0);
}
