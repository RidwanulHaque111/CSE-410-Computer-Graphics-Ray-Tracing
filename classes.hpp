#include <bits/stdc++.h>
using namespace std;
#include "vector.hpp"
int recursion_level;

//====================LIGHTS AND RAY CLASSES START==========================================
class Ray
{
public:
    Vector3D start;
    Vector3D dir;
    Ray(Vector3D, Vector3D);
};

Ray::Ray(Vector3D s, Vector3D d)
{
    start = s;
    dir = d;
    dir.normalize();
}

class PointLight
{
public:
    Vector3D light_pos;
    double color[3];
    double falloff;
    PointLight() {}

    virtual void print()
    {
        cout << "\nPoint Light==========\n";
        cout << light_pos;
        cout << "colors: " << color[0] << "," << color[1] << "," << color[2] << endl;
    }

    void read_pointlight(ifstream &ifs)
    {
        ifs >> light_pos >> falloff;
        // Create a random number generator engine
        std::random_device rd;  // Obtain a random seed from the hardware
        std::mt19937 gen(rd()); // Seed the generator

        // Create a distribution that produces random numbers between 0 and 1
        std::uniform_real_distribution<> dis(0.1, 1.000000);

        // Generate random colors between 0 to 1(inclusive)
        color[0] = dis(gen); //dis(gen);
        color[1] = dis(gen);
        color[2] = dis(gen);

    }

    virtual void draw()
    {
        glPushMatrix(); //==================
        glTranslated(light_pos.x, light_pos.y, light_pos.z);
        Vector3D points[100][100];
        int i,j;
        int slices = 30, stacks = 36;
        double h,r,length = 1.0;

        //generate points
        for(i=0; i<=stacks; i++)
        {
            h=length*sin(((double)i/(double)stacks)*(pi/2));
            r=length*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0; j<=slices; j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }

        // Calculate the distance from the light source to the current point
        double distance = sqrt(
                              pow(points[i][j].x - light_pos.x, 2) +
                              pow(points[i][j].y - light_pos.y, 2) +
                              pow(points[i][j].z - light_pos.z, 2)
                          );

        // Calculate the intensity based on falloff
        double intensity = 1.0 / pow(distance, falloff);
        //double intensity = 1.0;

        glColor3f(color[0] * intensity, color[1] * intensity, color[2] * intensity);



        //draw quads using generated points
        for(i=0; i<stacks; i++)
        {
            for(j=0; j<slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
                glEnd();
            }
        }

        glPopMatrix(); //===============
    }
};


class SpotLight: public PointLight
{
public:
    Vector3D light_dir;
    double cutoff_angle;

    SpotLight() {}
    void print() override
    {
        cout << "\nSpot Light==========\n";
        cout << light_pos;
        cout << "colors: " << color[0] << "," << color[1] << "," << color[2] << endl;
        cout << light_dir;
        cout << "cutoff angle: " << cutoff_angle << endl;
    }

    void read_spotlight(ifstream &ifs)
    {
        PointLight::read_pointlight(ifs);
        ifs >> light_dir >> cutoff_angle;
        light_dir.normalize();
    }

    void draw() override
    {
        PointLight::draw();
        Vector3D from = light_pos;
        Vector3D to = light_pos + light_dir*10.0;
        glBegin(GL_LINES);
        {
            glVertex3f(from.x, from.y, from.z);
            glVertex3f(to.x, to.y, to.z);
        }
        glEnd();
    }
};



//Global Light Data ========
vector <PointLight *> lights_list;

//====================LIGHTS AND RAY CLASSES END==========================================




//====================OBJECT CLASS START==========================================
class Object;
vector <Object *> objects;
double floor_size, tile_size;

class Object
{
public:
    double height, width, length;
    double color[3];
    double coEfficients[4];
    int shine;
    int obj_type;

    Object()
    {
        obj_type = 0;
    }
    virtual void draw() = 0;
    virtual Ray get_normal(Vector3D intersec_point, Ray incident_ray) = 0;
    void setShine(int);
    void print_obj();

    void read_obj(ifstream &);
    void setColor(double, double, double);
    void setCoEfficients(double, double, double, double);
    virtual double intersect(Ray ray, double *col, int level)
    {
        return -1.0;
    }


    // GENERIC PHONG LIGHTING FUNCTION FOR ALL OBJECTS =========================================
    void illuminate(Ray ray, double *col, double tMin, bool isFloor=false)
    {
        Vector3D intersec_point = ray.start + ray.dir * tMin;
        double floor_col = 0;
        int tile_row, tile_col;
        if(isFloor)
        {
            tile_row = (intersec_point.x + floor_size/2.0) / tile_size;
            tile_col = (intersec_point.y + floor_size/2.0) / tile_size;
            if((tile_col+tile_row)%2 == 1) floor_col = 1; //white if odd
            else return; //LIGHTING NOT DONE FOR DARK TILES!!! ==============================

            //Ambient Lighting for Floor
            for(int i=0; i<3; i++)
                col[i] += floor_col * coEfficients[0]; //AMBIENT
        }
        else
        {
            //Ambient Lighting for others
            for(int i=0; i<3; i++)
                col[i] += color[i] * coEfficients[0]; //AMBIENT
        }

        //Specular and Diffuse Lighting
        for (int i=0; i<lights_list.size(); i++)
        {
            bool light_obstructed = false;
            PointLight *light = lights_list.at(i);
            Vector3D incident_vect = intersec_point - light->light_pos;

            //if spotlight, check if object is within cutoff angle
            if(SpotLight* spotlight = dynamic_cast<SpotLight*>(light))
            {
                double angle = (incident_vect^spotlight->light_dir);
                angle /= (incident_vect.get_length()*spotlight->light_dir.get_length());
                angle = acos(angle);
                if(angle > rad(spotlight->cutoff_angle)) continue;
            }

            double light_distance = incident_vect.get_length();
            if (light_distance < ZERO) //light source inside object
                continue;

            //Ray L_ray = Ray of Incidence | N_ray = Ray of Normal | R_ray = Ray of Reflection
            Ray L_ray = Ray(light->light_pos, incident_vect);
            Ray N_ray = get_normal(intersec_point, L_ray);
            Vector3D R_dir = L_ray.dir - N_ray.dir * ((L_ray.dir^N_ray.dir)*2); //R = L - 2(L.N)N
            Ray R_ray = Ray(intersec_point, R_dir);

            //SHADOWS: is the incident ray being blocked by an object?
            for(int k=0; k<objects.size(); k++)
            {
                double t = objects.at(k)->intersect(L_ray, col, 0);
                if (t > 0 && t + ZERO < light_distance)
                {
                    light_obstructed = true;
                    break;
                }
            }

            if(light_obstructed) continue; //light
            for (int i = 0; i < 3; i++)
            {
                if(isFloor)
                {
                    //diffuse lighting component
                    col[i] += light->color[i] * coEfficients[1] * max((-L_ray.dir)^N_ray.dir, 0.0) * floor_col;
                    //specular lighting component
                    col[i] += light->color[i] * coEfficients[2] * pow(max((-ray.dir)^R_ray.dir, 0.0), shine);
                }
                else
                {
                    //diffuse lighting component
                    col[i] += light->color[i] * coEfficients[1] * max((-L_ray.dir)^N_ray.dir, 0.0) * color[i];
                    //specular lighting component
                    col[i] += light->color[i] * coEfficients[2] * pow(max((-ray.dir)^R_ray.dir, 0.0), shine);
                }
            }


        }// end of diffuse, specular loop ========================

        //keep colors within range:
        for(int i=0; i<3; i++) if(col[i] > 1.0) col[i] = 1.0;
        for(int i=0; i<3; i++) if(col[i] < ZERO) col[i] = 0.0;

    }// end of function ================================================================



    // GENERIC RECURSIVE REFLECTION FUNCTION FOR ALL OBJECTS =========================================
    void recursive_reflection(Ray ray, double *col, int level, double tMin)
    {
        //the last bounce of light is done
        if(level >= recursion_level) return;

        //Ray L_ray = ray; given as parameter | N_ray = Ray of Normal | R_ray = Ray of Reflection
        Vector3D intersec_point = ray.start + ray.dir * tMin;
        Ray N_ray = get_normal(intersec_point, ray);
        Vector3D R_dir = ray.dir - N_ray.dir * ((ray.dir^N_ray.dir)*2); //R = L - 2(L.N)N
        Ray R_ray = Ray(intersec_point, R_dir);
        R_ray.start = R_ray.start + R_ray.dir * ZERO;

        double nearest = INT_MAX;
        double tmin_temp = INT_MAX;
        for(int k=0; k<objects.size(); k++)
        {
            double t = objects.at(k)->intersect(R_ray, col, 0);
            if(t>0.0 && t<tmin_temp)
            {
                tmin_temp = t;
                nearest = k;
            }
        }

        if (nearest == INT_MAX) return; //no object reflects the light :(

        //empty color array to find reflected color
        double *color_temp = new double[3];
        for (int i = 0; i < 3; i++) color_temp[i] = 0;

        //sets the reflected color
        tmin_temp = objects[nearest]->intersect(R_ray, color_temp, level + 1);
        for (int i = 0; i < 3; i++) col[i] += color_temp[i] * coEfficients[3];

        //always clear memory :-D
        delete[] color_temp;
    }

};

// ============================================


//Global camera data =========
Vector3D c_pos, u, l, r;
// ================================


void Object::setColor(double r, double g, double b)
{
    color[0] = r;
    color[1] = g;
    color[2] = b;
}

void Object::setShine(int s)
{
    shine = s;
}

void Object::setCoEfficients(double amb, double diff, double spec, double refl)
{
    coEfficients[0] = amb;
    coEfficients[1] = diff;
    coEfficients[2] = spec;
    coEfficients[3] = refl;
}

void Object::read_obj(ifstream &ifs)
{
    ifs >> color[0] >> color[1] >> color[2];
    ifs >> coEfficients[0] >> coEfficients[1] >> coEfficients[2] >> coEfficients[3];
    ifs >> shine;
}

void Object::print_obj()
{
    cout << "colors: " << color[0] << "," << color[1] << "," << color[2] << endl;
    cout << "coeff: " << coEfficients[0] << "," << coEfficients[1] << ",";
    cout << coEfficients[2] << "," << coEfficients[3] << endl;
    cout << "shine: " << shine << endl;
}
//====================OBJECT CLASS END==========================================


//====================SPHERE CLASS START==========================================
class Sphere : public Object
{
public:
    Vector3D centre;

    Sphere()
    {
        obj_type = 1;
    }
    void print();
    void read_sphere(ifstream &);

    virtual void draw()
    {
        glPushMatrix(); //=========================================
        glTranslated(centre.x, centre.y, centre.z);
        Vector3D points[100][100];
        int i,j;
        int slices = 30, stacks = 32;
        double h,r;
        glColor3f(color[0], color[1], color[2]);

        //generate points
        for(i=0; i<=stacks; i++)
        {
            h=length*sin(((double)i/(double)stacks)*(pi/2));
            r=length*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0; j<=slices; j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0; i<stacks; i++)
        {
            //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
            for(j=0; j<slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
                glEnd();
            }
        }

        glPopMatrix(); //=====================================================
    }

    double get_tMin(Ray ray)
    {
        ray.start = ray.start - centre; //translate using centre
        double td = -ray.start^ray.dir; //sign of cosine of Ro and Rd
        if(td < 0.0) return -1.0;

        double d_squared = (ray.start^ray.start) - td*td;
        if(d_squared > length*length) return -1.0; //length is radius

        double offset = sqrt(length*length - d_squared);
        double t = td - offset;
        if(t > 0.0) return t;

        t = td + offset;
        if(t > 0.0) return t;

        return -1.0;
    }

    virtual Ray get_normal(Vector3D intersec_point, Ray incident_ray)
    {
        return Ray(intersec_point, intersec_point - centre);
    }

    virtual double intersect(Ray ray, double *col, int level)
    {
        double tMin = get_tMin(ray);
        if(level == 0) return tMin;

        //illumination function same for all
        illuminate(ray, col, tMin);

        //recursive reflection of light
        recursive_reflection(ray, col, level, tMin);

        //return min val of intersection equation parameter: t
        return tMin;
    }

};

void Sphere::read_sphere(ifstream &ifs)
{
    ifs >> centre >> length;
    read_obj(ifs);
}

void Sphere::print()
{
    cout << "\nSphere=================\n";
    cout << "centre: " << centre;
    cout << "rad: " << length << endl;
    print_obj();
}
//====================SPHERE CLASS END==========================================


//=======================CUBE CLASS START==========================================
class Cube : public Object
{
public:
    Vector3D lowerLeft; // Bottom lower left point of cube
    double side; // Cube side length

    Cube()
    {
        obj_type = 2;
    }

    void print();
    void read_cube(ifstream&);

    virtual void draw()
    {
        glPushMatrix();
        glTranslated(lowerLeft.x + side / 2, lowerLeft.y + side / 2, lowerLeft.z + side / 2);
        glColor3f(color[0], color[1], color[2]);

        // Draw the front face
        glBegin(GL_QUADS);
        glVertex3d(-side / 2, -side / 2, side / 2);
        glVertex3d(side / 2, -side / 2, side / 2);
        glVertex3d(side / 2, side / 2, side / 2);
        glVertex3d(-side / 2, side / 2, side / 2);
        glEnd();

        // Draw the back face
        glBegin(GL_QUADS);
        glVertex3d(side / 2, -side / 2, -side / 2);
        glVertex3d(-side / 2, -side / 2, -side / 2);
        glVertex3d(-side / 2, side / 2, -side / 2);
        glVertex3d(side / 2, side / 2, -side / 2);
        glEnd();

        // Draw the left face
        glBegin(GL_QUADS);
        glVertex3d(-side / 2, -side / 2, -side / 2);
        glVertex3d(-side / 2, -side / 2, side / 2);
        glVertex3d(-side / 2, side / 2, side / 2);
        glVertex3d(-side / 2, side / 2, -side / 2);
        glEnd();

        // Draw the right face
        glBegin(GL_QUADS);
        glVertex3d(side / 2, -side / 2, side / 2);
        glVertex3d(side / 2, -side / 2, -side / 2);
        glVertex3d(side / 2, side / 2, -side / 2);
        glVertex3d(side / 2, side / 2, side / 2);
        glEnd();

        // Draw the top face
        glBegin(GL_QUADS);
        glVertex3d(-side / 2, side / 2, side / 2);
        glVertex3d(side / 2, side / 2, side / 2);
        glVertex3d(side / 2, side / 2, -side / 2);
        glVertex3d(-side / 2, side / 2, -side / 2);
        glEnd();

        // Draw the bottom face
        glBegin(GL_QUADS);
        glVertex3d(-side / 2, -side / 2, -side / 2);
        glVertex3d(side / 2, -side / 2, -side / 2);
        glVertex3d(side / 2, -side / 2, side / 2);
        glVertex3d(-side / 2, -side / 2, side / 2);
        glEnd();

        glPopMatrix();
    }


    double get_tMin(Ray ray)
    {
        Vector3D dirfrac(1.0 / ray.dir.x, 1.0 / ray.dir.y, 1.0 / ray.dir.z);
        double t1 = (lowerLeft.x - ray.start.x) * dirfrac.x;
        double t2 = (lowerLeft.x + side - ray.start.x) * dirfrac.x;
        double t3 = (lowerLeft.y - ray.start.y) * dirfrac.y;
        double t4 = (lowerLeft.y + side - ray.start.y) * dirfrac.y;
        double t5 = (lowerLeft.z - ray.start.z) * dirfrac.z;
        double t6 = (lowerLeft.z + side - ray.start.z) * dirfrac.z;

        double tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
        double tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

        // If tmax < 0, ray is intersecting AABB, but the whole AABB is behind us
        if (tmax < 0)
        {
            return -1.0;
        }

        // If tmin > tmax, ray doesn't intersect AABB
        if (tmin > tmax)
        {
            return -1.0;
        }

        // Else, ray intersects AABB
        return tmin;
    }

    virtual Ray get_normal(Vector3D intersec_point, Ray incident_ray)
    {
        Vector3D center = lowerLeft + Vector3D(side / 2, side / 2, side / 2);

        // Determine which face was hit based on the intersection point
        Vector3D diff = intersec_point - center;
        double absX = fabs(diff.x);
        double absY = fabs(diff.y);
        double absZ = fabs(diff.z);

        Vector3D normal;
        if (absX > absY && absX > absZ)
        {
            normal = Vector3D(diff.x, 0, 0);
        }
        else if (absY > absX && absY > absZ)
        {
            normal = Vector3D(0, diff.y, 0);
        }
        else
        {
            normal = Vector3D(0, 0, diff.z);
        }

        // Make sure the normal vector points outwards
        if (-(normal^incident_ray.dir) >= 0)
        {
            normal = -normal;
        }

        return Ray(intersec_point, normal);
    }

    virtual double intersect(Ray ray, double *col, int level)
    {
        double tMin = get_tMin(ray);
        if (level == 0) return tMin;

        // Illuminate function same for all
        illuminate(ray, col, tMin);

        // Recursive reflection of light
        recursive_reflection(ray, col, level, tMin);

        // Return min value of intersection equation parameter: t
        return tMin;
    }
};

void Cube::read_cube(ifstream &ifs)
{
    ifs >> lowerLeft >> side;
    read_obj(ifs);
}

void Cube::print()
{
    cout << "\nCube===================\n";
    cout << "Lower Left: " << lowerLeft;
    cout << "Side: " << side << endl;
    print_obj();
}
//=======================CUBE CLASS END==========================================




//=======================PYRAMID CLASS START==========================================
class Pyramid : public Object
{
public:
    Vector3D lowestPoint;// Lowest point of a quad pyramid which is the intersection of the diagonals of the base
    double width, height; // Pyramid width and height

    Pyramid()
    {
        obj_type = 3;
    }

    void print();
    void read_pyramid(ifstream&);

    virtual void draw()
    {
        glPushMatrix();

        // Calculate the apex of the pyramid
        Vector3D apex = lowestPoint + Vector3D(0, 0, height);

        // Draw the base of the pyramid
        glBegin(GL_QUADS);
        glColor3f(color[0], color[1], color[2]);
        glVertex3d(lowestPoint.x, lowestPoint.y, lowestPoint.z);
        glVertex3d(lowestPoint.x + width, lowestPoint.y, lowestPoint.z);
        glVertex3d(lowestPoint.x + width, lowestPoint.y + width, lowestPoint.z);
        glVertex3d(lowestPoint.x, lowestPoint.y + width, lowestPoint.z);
        glEnd();

        // Draw the triangular faces of the pyramid
        glBegin(GL_TRIANGLES);
        glColor3f(color[0], color[1], color[2]);
        // First triangular face
        glVertex3d(lowestPoint.x + width / 2, lowestPoint.y + width / 2, lowestPoint.z + height);
        glVertex3d(lowestPoint.x, lowestPoint.y, lowestPoint.z);
        glVertex3d(lowestPoint.x + width, lowestPoint.y, lowestPoint.z);
        // Second triangular face
        glVertex3d(lowestPoint.x + width / 2, lowestPoint.y + width / 2, lowestPoint.z + height);
        glVertex3d(lowestPoint.x + width, lowestPoint.y, lowestPoint.z);
        glVertex3d(lowestPoint.x + width, lowestPoint.y + width, lowestPoint.z);
        // Third triangular face
        glVertex3d(lowestPoint.x + width / 2, lowestPoint.y + width / 2, lowestPoint.z + height);
        glVertex3d(lowestPoint.x + width, lowestPoint.y + width, lowestPoint.z);
        glVertex3d(lowestPoint.x, lowestPoint.y + width, lowestPoint.z);
        // Fourth triangular face
        glVertex3d(lowestPoint.x + width / 2, lowestPoint.y + width / 2, lowestPoint.z + height);
        glVertex3d(lowestPoint.x, lowestPoint.y + width, lowestPoint.z);
        glVertex3d(lowestPoint.x, lowestPoint.y, lowestPoint.z);
        glEnd();

        glPopMatrix();
    }

    double get_tMin(Ray ray)
    {
        // Calculate the tMin value for intersection with the pyramid
        Vector3D normal = unit_vector(cross_product(Vector3D(width, 0, 0), Vector3D(0, width, 0)));
        Vector3D v0 = lowestPoint;
        Vector3D v1 = v0 + Vector3D(width, 0, 0);
        Vector3D v2 = v1 + Vector3D(0, width, 0);
        Vector3D v3 = v0 + Vector3D(0, width, 0);
        Vector3D v4 = v0 + Vector3D(0, 0, height);

        double tMin = INT_MAX;

        // Check intersection with the base of the pyramid
        double tBase = (v0 - ray.start) ^ normal / (ray.dir ^ normal);
        if (tBase > ZERO)
        {
            Vector3D intersection_point = ray.start + ray.dir * tBase;
            if (intersection_point.x >= v0.x && intersection_point.x <= v1.x &&
                    intersection_point.y >= v0.y && intersection_point.y <= v3.y)
            {
                tMin = min(tMin, tBase);
            }
        }

        // Check intersection with the four triangular faces of the pyramid
        Vector3D points[4] = {v0, v1, v2, v3};
        for (int i = 0; i < 4; i++)
        {
            Vector3D v1 = points[i];
            Vector3D v2 = v4;
            Vector3D v3 = points[(i + 1) % 4];

            Vector3D edge1 = v2 - v1;
            Vector3D edge2 = v3 - v1;
            Vector3D h = ray.dir.cross(edge2);
            double a = edge1 ^ h;

            if (a > -ZERO && a < ZERO)
                continue;

            double f = 1.0 / a;
            Vector3D s = ray.start - v1;
            double u = f * (s ^ h);

            if (u < 0.0 || u > 1.0)
                continue;

            Vector3D q = s.cross(edge1);
            double v = f * (ray.dir ^ q);

            if (v < 0.0 || u + v > 1.0)
                continue;

            double t = f * (edge2 ^ q);

            if (t > ZERO)
            {
                tMin = min(tMin, t);
            }
        }

        return tMin;
    }

    virtual Ray get_normal(Vector3D intersec_point, Ray incident_ray)
    {
        // Calculate the normal vector at the intersection point on the pyramid
        Vector3D normal = unit_vector(cross_product(Vector3D(width, 0, 0), Vector3D(0, width, 0)));
        return Ray(intersec_point, normal);

    }


    virtual double intersect(Ray ray, double *col, int level)
    {
        double tMin = get_tMin(ray);
        if (level == 0)
            return tMin;

        // Check if there was an intersection with the pyramid
        if (tMin < INFINITY)
        {
            // Calculate the intersection point
            Vector3D intersection_point = ray.start + ray.dir * tMin;

            // Calculate the normal vector at the intersection point
            Vector3D normal;

            // Calculate the normal based on the face of the pyramid that was hit
            Vector3D apex = lowestPoint + Vector3D(width / 2, width / 2, height);
            Vector3D base1 = lowestPoint;
            Vector3D base2 = lowestPoint + Vector3D(width, 0, 0);

            if (intersection_point.z > apex.z - ZERO)
            {
                // The intersection point is on the base of the pyramid
                normal = unit_vector(cross_product(base2 - base1, Vector3D(0, 0, -1)));
            }
            else
            {
                // The intersection point is on one of the triangular faces
                Vector3D base_mid = lowestPoint + Vector3D(width / 2, width / 2, 0);

                if (intersection_point == apex)
                {
                    // Special case: intersection point is at the apex
                    normal = unit_vector(apex - base_mid);
                }
                else
                {
                    // Calculate the normal by taking the cross product of two edges of the triangular face
                    Vector3D edge1, edge2;
                    if (intersection_point == base1)
                    {
                        edge1 = base2 - apex;
                        edge2 = intersection_point - apex;
                    }
                    else if (intersection_point == base2)
                    {
                        edge1 = base1 - apex;
                        edge2 = intersection_point - apex;
                    }
                    else
                    {
                        edge1 = base1 - apex;
                        edge2 = base2 - apex;
                    }
                    normal = unit_vector(cross_product(edge1, edge2));
                }
            }

            // Calculate the reflected ray direction
            Vector3D reflected_dir = ray.dir - (ray.dir * normal) * normal * 2.0;

            // Set the color of the pyramid at the intersection point
            col[0] = color[0] * coEfficients[0];
            col[1] = color[1] * coEfficients[0];
            col[2] = color[2] * coEfficients[0];
            

            // Recursive reflection
            Ray reflected_ray(intersection_point, reflected_dir);
            double reflected_color[3] = {0.0, 0.0, 0.0};
            double tReflect = intersect(reflected_ray, reflected_color, level - 1);

            // Add the reflected color to the current color
            col[0] += 0.2 * reflected_color[0];
            col[1] += 0.2 * reflected_color[1];
            col[2] += 0.2 * reflected_color[2];

            return tMin;
        }
        return -1.0; // No intersection
    }


};

void Pyramid::read_pyramid(ifstream &ifs)
{
    ifs >> lowestPoint >> width >> height;
    read_obj(ifs);
}

void Pyramid::print()
{
    cout << "\nPyramid===================\n";
    cout << "Lowest Point: " << lowestPoint;
    cout << "Width: " << width << endl;
    cout << "Height: " << height << endl;
    print_obj();
}
//=======================PYRAMID CLASS END==========================================




//=======================FLOOR CLASS START==========================================
class Floor : public Object
{
public:
    Floor(int floorSize, int tileSize)
    {
        obj_type = 4;
        length = floorSize;
        width = tileSize;
        //global
        floor_size = floorSize;
        tile_size = tileSize;
    }

    virtual void draw()
    {
        int row = length / width, col = length / width;
        int quad_length = length / 2;

        int shiftX = 0; // Shift amount for X axis
        int shiftY = 0; // Shift amount for Y axis

        // Calculate shift amounts based on camera movement
        if (c_pos.x < 0)
        {
            shiftX = static_cast<int>(abs(c_pos.x) / width);
        }
        if (c_pos.y < 0)
        {
            shiftY = static_cast<int>(abs(c_pos.y) / width);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                // xy-plane so z = 0
                Vector3D tile_start_point(width * (j - shiftX) - quad_length, width * (i - shiftY) - quad_length, 0);

                // setting the color
                if ((i + j + shiftX + shiftY) % 2 == 0)
                    glColor3f(0, 0, 0);
                else
                    glColor3f(color[0], color[1], color[2]);

                glBegin(GL_QUADS);
                {
                    glVertex3f(tile_start_point.x, tile_start_point.y, 0);
                    glVertex3f(tile_start_point.x + width, tile_start_point.y, 0);
                    glVertex3f(tile_start_point.x + width, tile_start_point.y + width, 0);
                    glVertex3f(tile_start_point.x, tile_start_point.y + width, 0);
                }
                glEnd();
            }
        }
    }

    virtual Ray get_normal(Vector3D intersec_point, Ray incident_ray)
    {
        Vector3D normal(0, 0, 1);
        if (incident_ray.start.z > 0) return Ray(intersec_point, normal);
        return Ray(intersec_point, -normal);
    }

    double get_tMin(Ray ray)
    {
        Vector3D n(0, 0, 1);
        if(c_pos.z < 0.0) n = -n;

        double denom = n^ray.dir;
        if(denom == 0.0) return -1.0;

        //Po is orgin so ZERO!!
        double tMin = (-n^ray.start)/(denom);
        Vector3D intersec_point = ray.start + ray.dir * tMin;
        if (fabs(intersec_point.y) > length/2.0) return -1.0;   //length is floorSize :-D
        if (fabs(intersec_point.x) > length/2.0) return -1.0;
        return tMin;
    }


    virtual double intersect(Ray ray, double *col, int level)
    {
        double tMin = get_tMin(ray);
        if(level == 0) return tMin;

        //illumination function same for all
        illuminate(ray, col, tMin, true);

        //recursive reflection of light
        recursive_reflection(ray, col, level, tMin);

        //return min val of intersection equation parameter: t
        return tMin;
    }

};
//=======================FLOOR CLASS END==========================================
