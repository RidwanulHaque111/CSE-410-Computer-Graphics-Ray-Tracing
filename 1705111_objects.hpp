#include <bits/stdc++.h>
using namespace std;
#include "1705111_ray.hpp"

//global var
VectorPoint camera_pos, u, l, r;
int recursion_level;
int lightCounter = 1;

class PointLightObj
{
public:
    VectorPoint light_pos;
    VectorPoint temp_pos;
    double color[3];
    double falloff;

    PointLightObj() {}

    virtual void print()
    {
        cout << "\Light "<<lightCounter<<" Loaded."<<endl;
        lightCounter++;

    }

    void input_pointlight(ifstream &ifs)
    {
        ifs >> light_pos >> falloff;
        print();

        double falloffUpdate = falloff;
        update(falloffUpdate);

        // Create a random number generator engine
        std::random_device rd;  // Obtain a random seed from the hardware
        std::mt19937 gen(rd()); // Seed the generator

        // Create a distribution that produces random numbers between 0 and 1
        std::uniform_real_distribution<> dis(0.1, 1);

        // Generate random colors between 0 to 1(inclusive)
        color[0] = dis(gen); //dis(gen);
        color[1] = dis(gen);
        color[2] = dis(gen);
    }

    void setPosition(const VectorPoint& position)
    {
        light_pos = position;
    }

    void setColors(double red, double green, double blue)
    {
        color[0] = red;
        color[1] = green;
        color[2] = blue;
    }

    void setFalloff(double value)
    {
        falloff = value;
    }

    VectorPoint getPosition() const
    {
        return light_pos;
    }

    double* getColor()
    {
        return color;
    }

    double getFalloff() const
    {
        return falloff;
    }

    void update(double val)
    {
        VectorPoint updatedPoint = light_pos;
        updatedPoint.x += 0.1 * val;
        updatedPoint.y += 0.1 * val;
    }


    virtual void draw()
    {
        glPushMatrix(); //==================
        glTranslated(light_pos.x, light_pos.y, light_pos.z);
        VectorPoint points[100][100];
        int i,j;
        int slices = 30, stacks = 36;
        double h,r,length = 1.0;

        //generate points
        for(i=0; i<=stacks; i++)
        {
            r=length*cos(((double)i/(double)stacks)*(pi/2));
            h=length*sin(((double)i/(double)stacks)*(pi/2));

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


class SpotLight: public PointLightObj
{
public:
    double cutoff_angle;
    VectorPoint light_dir;
    VectorPoint temp_dir;

    SpotLight() {}
    void print() override
    {
        //cout<<"(Spot Light)"<<endl;
    }

    void input_spotlight(ifstream &ifs)
    {
        PointLightObj::input_pointlight(ifs);
        ifs >> light_dir >> cutoff_angle;
        light_dir.normalize();
        print();
    }

    void setLightDirection(const VectorPoint& direction)
    {
        light_dir = direction;
        light_dir.normalize();
    }

    VectorPoint getLightDirection() const
    {
        return light_dir;
    }

    void setCutoffAngle(double angle)
    {
        cutoff_angle = angle;
    }

    double getCutoffAngle() const
    {
        return cutoff_angle;
    }

    void setLightPosition(const VectorPoint& position)
    {
        light_pos = position;
    }

    void adjustIntensity(double factor)
    {
        color[0] *= factor;
        color[1] *= factor;
        color[2] *= factor;
    }


    void resetToDefault()
    {
        adjustIntensity(1.0);
        print();
        //light_pos = defaultLightPos;
        //color = defaultColor;
        //light_dir = defaultDirection;
        // cutoff_angle = defaultCutoffAngle;

    }


    void draw() override
    {

        PointLightObj::draw();
        VectorPoint from = light_pos;
        VectorPoint to = light_pos + light_dir*10.0;

        glBegin(GL_LINES);
        {
            glVertex3f(from.x, from.y, from.z);
            glVertex3f(to.x, to.y, to.z);
        }
        glEnd();
    }
};


vector <PointLightObj *> lights_list;


class ObjectClass;
vector <ObjectClass *> objects;
double floor_count, tile_count, floor_size, tile_size = 0.0;

class ObjectClass
{
public:
    double color[3];
    double coEfficients[4];
    int shine;
    int obj_type;
    double height, width, length;

    ObjectClass()
    {
        obj_type = 0;
    }

    void setShine(int);
    void print_object_details();
    void input_obj(ifstream &);
    virtual void draw() = 0;
    virtual Ray find_normal(VectorPoint intersecPoint, Ray incident_ray) = 0;
    void setColors(double, double, double);
    void setCoEfficients(double, double, double, double);

    virtual double find_intersect(Ray ray, double *col, int level)
    {
        return -1.0; //default
    }


    // GENERIC PHONG LIGHTING FUNCTION FOR ALL OBJECTS =========================================
    void illuminate(Ray ray, double *col, double tMin, bool isFloor=false)
    {
        int floorCounter = 0;
        VectorPoint intersecPoint = ray.start + ray.dir * tMin;
        double floor_col = 0;
        int tile_row, tile_col;
        if(isFloor)
        {
            tile_row = (intersecPoint.x + floor_size/2.0) / tile_size;
            tile_col = (intersecPoint.y + floor_size/2.0) / tile_size;
            if((tile_col+tile_row)%2 == 1) floor_col = 1; //white if odd
            else return;

            //Ambient Lighting for Floor
            for(int i=0; i<3; i++)
                col[i] += floor_col * coEfficients[0];
        }
        else
        {
            //Ambient Lighting for others
            for(int i=0; i<3; i++)
                col[i] += color[i] * coEfficients[0];
        }

        //Specular and Diffuse Lighting
        for (int i=0; i<lights_list.size(); i++)
        {
            bool isLightObstructed = false;
            PointLightObj *light = lights_list.at(i);
            VectorPoint incident_vect = intersecPoint - light->light_pos;

            //if spotlight, check if ObjectClass is within cutoff angle
            if(SpotLight* spotlight = dynamic_cast<SpotLight*>(light))
            {
                double angle = (incident_vect^spotlight->light_dir);
                angle /= (incident_vect.get_length()*spotlight->light_dir.get_length());
                angle = acos(angle);

                double currentAngle = angle;
                if(currentAngle > rad(spotlight->cutoff_angle)) continue;
            }

            double light_distance = incident_vect.get_length();
            if (light_distance < ZERO) //light source inside ObjectClass
                continue;

            //Ray L_ray = Ray of Incidence | N_ray = Ray of Normal | R_ray = Ray of Reflection
            Ray L_ray = Ray(light->light_pos, incident_vect);
            Ray N_ray = find_normal(intersecPoint, L_ray);
            VectorPoint R_dir = L_ray.dir - N_ray.dir * ((L_ray.dir^N_ray.dir)*2); //R = L - 2(L.N)N
            Ray R_ray = Ray(intersecPoint, R_dir);

            //SHADOWS: is the incident ray being blocked by an ObjectClass?
            for(int k=0; k<objects.size(); k++)
            {
                double t = objects.at(k)->find_intersect(L_ray, col, 0);
                if (t > 0 && t + ZERO < light_distance)
                {
                    isLightObstructed = true;
                    break;
                }
            }

            if(isLightObstructed) continue; //light
            for (int i = 0; i < 3; i++)
            {
                if(isFloor)
                {
                    floorCounter = floorCounter + 1;
                    //diffuse lighting components
                    col[i] += light->color[i] * coEfficients[1] * max((-L_ray.dir)^N_ray.dir, 0.0) * floor_col;
                    //specular lighting components
                    col[i] += light->color[i] * coEfficients[2] * pow(max((-ray.dir)^R_ray.dir, 0.0), shine);
                }
                else
                {
                    //diffuse lighting components
                    col[i] += light->color[i] * coEfficients[1] * max((-L_ray.dir)^N_ray.dir, 0.0) * color[i];
                    //specular lighting components
                    col[i] += light->color[i] * coEfficients[2] * pow(max((-ray.dir)^R_ray.dir, 0.0), shine);
                }
            }

        }

        //keep colors within range:
        for(int i=0; i<3; i++) if(col[i] > 1.0) col[i] = 1.0;
        for(int i=0; i<3; i++) if(col[i] < ZERO) col[i] = 0.0;

    }// end of function


    // GENERIC RECURSIVE REFLECTION FUNCTION FOR ALL OBJECTS ===================
    void rec_reflection(Ray ray, double *col, int level, double tMin)
    {
        //the last bounce of light is done
        if(level >= recursion_level) return;

        //Ray L_ray = ray; given as parameter | N_ray = Ray of Normal | R_ray = Ray of Reflection
        VectorPoint intersecPoint = ray.start + ray.dir * tMin;
        Ray N_ray = find_normal(intersecPoint, ray);
        VectorPoint R_dir = ray.dir - N_ray.dir * ((ray.dir^N_ray.dir)*2); //R = L - 2(L.N)N
        Ray R_ray = Ray(intersecPoint, R_dir);
        R_ray.start = R_ray.start + R_ray.dir * ZERO;

        double tmin_temp = INT_MAX;
        double nearest = INT_MAX;

        for(int k=0; k<objects.size(); k++)
        {
            double t = objects.at(k)->find_intersect(R_ray, col, 0);
            if(t>0.0 && t<tmin_temp)
            {
                tmin_temp = t;
                nearest = k;
            }
        }

        if (nearest == INT_MAX) return; //no ObjectClass reflects the lightts

        //empty color array to find reflected color
        double *color_temp = new double[3];
        for (int i = 0; i < 3; i++) color_temp[i] = 0;

        //setting the reflected color
        tmin_temp = objects[nearest]->find_intersect(R_ray, color_temp, level + 1);
        for (int i = 0; i < 3; i++) col[i] += color_temp[i] * coEfficients[3];


        //clear memory
        delete[] color_temp;
    }

};


void ObjectClass::setColors(double r, double g, double b)
{
    color[0] = r;
    color[1] = g;
    color[2] = b;
}

void ObjectClass::setShine(int s)
{
    shine = s;
}

void ObjectClass::setCoEfficients(double amb, double diff, double spec, double refl)
{
    coEfficients[0] = amb;
    coEfficients[1] = diff;
    coEfficients[2] = spec;
    coEfficients[3] = refl;
}

void ObjectClass::input_obj(ifstream &ifs)
{
    ifs >> color[0] >> color[1] >> color[2];
    ifs >> coEfficients[0] >> coEfficients[1] >> coEfficients[2] >> coEfficients[3];
    ifs >> shine;
}

void ObjectClass::print_object_details()
{
    cout << "colors: " << color[0] << "," << color[1] << "," << color[2] << endl;

}




class Floor : public ObjectClass
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
        if (camera_pos.x < 0)
        {
            shiftX = static_cast<int>(abs(camera_pos.x) / width);
        }
        if (camera_pos.y < 0)
        {
            shiftY = static_cast<int>(abs(camera_pos.y) / width);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                // xy-plane so z = 0
                VectorPoint tile_start_point(width * (j - shiftX) - quad_length, width * (i - shiftY) - quad_length, 0);

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

    virtual Ray find_normal(VectorPoint intersecPoint, Ray incident_ray)
    {
        VectorPoint normal(0, 0, 1);
        if (incident_ray.start.z > 0) return Ray(intersecPoint, normal);
        return Ray(intersecPoint, -normal);
    }

    double find_tMin(Ray ray)
    {
        VectorPoint n(0, 0, 1);
        if(camera_pos.z < 0.0) n = -n;

        double denom = n^ray.dir;
        if(denom == 0.0) return -1.0;

        //Po is orgin so ZERO!!
        double tMin = (-n^ray.start)/(denom);
        VectorPoint intersecPoint = ray.start + ray.dir * tMin;
        if (fabs(intersecPoint.y) > length/2.0) return -1.0;
        if (fabs(intersecPoint.x) > length/2.0) return -1.0;
        return tMin;
    }

    virtual double find_intersect(Ray ray, double *col, int level)
    {
        double tMin = find_tMin(ray);
        if(level == 0) return tMin;

        //illumination function same for all
        illuminate(ray, col, tMin, true);

        //recursive reflection of light
        rec_reflection(ray, col, level, tMin);

        //return min val of find_intersection equation parameter: t
        return tMin;
    }

};




class Sphere : public ObjectClass
{
public:
    VectorPoint centre;
    double Radius;

    Sphere()
    {
        obj_type = 1;
    }

    void print();
    void input_sphere(ifstream &);


    void setCenter(const VectorPoint& newCenter)
    {
        centre = newCenter;
    }

    // Function to set the sphere's radius
    void setRadius(double newRadius)
    {
        length = newRadius;
    }

    // Function to get the sphere's center
    VectorPoint getCenter() const
    {
        return centre;
    }

    // Function to get the sphere's radius
    double getRadius() const
    {
        return length;
    }

    // Function to calculate the sphere's volume
    double calculateVolume() const
    {
        // Assuming the sphere is made of a uniform material
        return (4.0 / 3.0) * 3.14159265359 * length * length * length;
    }

    virtual void draw()
    {
        glPushMatrix();
        glTranslated(centre.x, centre.y, centre.z);
        VectorPoint CENTRE = getCenter();
        VectorPoint points[100][100];
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

        glPopMatrix();
    }

    double find_tMin(Ray ray)
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

    virtual Ray find_normal(VectorPoint intersecPoint, Ray incident_ray)
    {
        VectorPoint IntersectionPoint = intersecPoint;
        return Ray(IntersectionPoint, IntersectionPoint - centre);
    }

    virtual double find_intersect(Ray ray, double *col, int level)
    {
        double tMin = find_tMin(ray);
        if(level == 0) return tMin;

        //illumination function same for all
        illuminate(ray, col, tMin);

        //recursive reflection of light
        rec_reflection(ray, col, level, tMin);

        //return min val of find_intersection equation parameter: t
        return tMin;
    }
};

void Sphere::input_sphere(ifstream &ifs)
{
    ifs >> centre >> length;
    input_obj(ifs);
}

void Sphere::print()
{
    print_object_details();
}


class Cube : public ObjectClass
{
public:
    VectorPoint lowerLeft; // Bottom lower left point of cube
    double side; // Cube side length

    Cube()
    {
        obj_type = 2;
    }

    void print();
    void input_cube(ifstream&);


    double calculateVolume()
    {
        return side * side * side;
    }

    // New function 2: Calculate and return the surface area of the cube
    double calculateSurfaceArea()
    {
        return 6.0 * side * side;
    }

    // New function 3: Scale the cube by a given factor
    void scale(double scaleFactor)
    {
        side *= scaleFactor;
    }


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


    double find_tMin(Ray ray)
    {
        VectorPoint dirfrac(1.0 / ray.dir.x, 1.0 / ray.dir.y, 1.0 / ray.dir.z);
        double t1 = (lowerLeft.x - ray.start.x) * dirfrac.x;
        double t2 = (lowerLeft.x + side - ray.start.x) * dirfrac.x;
        double t3 = (lowerLeft.y - ray.start.y) * dirfrac.y;
        double t4 = (lowerLeft.y + side - ray.start.y) * dirfrac.y;
        double t5 = (lowerLeft.z - ray.start.z) * dirfrac.z;
        double t6 = (lowerLeft.z + side - ray.start.z) * dirfrac.z;

        double tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
        double tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

        // If tmax < 0, ray is find_intersecting AABB, but the whole AABB is behind us
        if (tmax < 0)
        {
            return -1.0;
        }

        // If tmin > tmax, ray doesn't find_intersect AABB
        if (tmin > tmax)
        {
            return -1.0;
        }

        // Else, ray find_intersects AABB
        return tmin;
    }

    virtual Ray find_normal(VectorPoint intersecPoint, Ray incident_ray)
    {
        VectorPoint center = lowerLeft + VectorPoint(side / 2, side / 2, side / 2);

        // Determine which face was hit based on the find_intersection point
        VectorPoint diff = intersecPoint - center;
        double absX = fabs(diff.x);
        double absY = fabs(diff.y);
        double absZ = fabs(diff.z);

        VectorPoint normal;
        if (absX > absY && absX > absZ)
        {
            normal = VectorPoint(diff.x, 0, 0);
        }
        else if (absY > absX && absY > absZ)
        {
            normal = VectorPoint(0, diff.y, 0);
        }
        else
        {
            normal = VectorPoint(0, 0, diff.z);
        }

        // Make sure the normal vector points outwards
        if (-(normal^incident_ray.dir) >= 0)
        {
            normal = -normal;
        }

        return Ray(intersecPoint, normal);
    }

    virtual double find_intersect(Ray ray, double *col, int level)
    {
        double tMin = find_tMin(ray);
        if (level == 0) return tMin;

        // Illuminate function same for all
        illuminate(ray, col, tMin);

        // Recursive reflection of light
        rec_reflection(ray, col, level, tMin);

        // Return min value of find_intersection equation parameter: t
        return tMin;
    }


};

void Cube::input_cube(ifstream &ifs)
{
    ifs >> lowerLeft >> side;
    input_obj(ifs);
}

void Cube::print()
{
    print_object_details();
}


class Pyramid : public ObjectClass
{
public:
    VectorPoint lowestPoint;// Lowest point of a quad pyramid which is the find_intersection of the diagonals of the base
    double width, height; // Pyramid width and height

    Pyramid()
    {
        obj_type = 3;
    }

    void print();
    void input_pyramid(ifstream&);

    virtual void draw()
    {
        glPushMatrix();

        // Calculate the apex of the pyramid
        VectorPoint apex = lowestPoint + VectorPoint(0, 0, height);

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

    double find_tMin(Ray ray)
    {
        // Calculate the tMin value for find_intersection with the pyramid
        VectorPoint normal = unit_vector(cross_product(VectorPoint(width, 0, 0), VectorPoint(0, width, 0)));
        VectorPoint v0 = lowestPoint;
        VectorPoint v1 = v0 + VectorPoint(width, 0, 0);
        VectorPoint v2 = v1 + VectorPoint(0, width, 0);
        VectorPoint v3 = v0 + VectorPoint(0, width, 0);
        VectorPoint v4 = v0 + VectorPoint(0, 0, height);

        double tMin = INT_MAX;

        // Check find_intersection with the base of the pyramid
        double tBase = (v0 - ray.start) ^ normal / (ray.dir ^ normal);
        if (tBase > ZERO)
        {
            VectorPoint find_intersection_point = ray.start + ray.dir * tBase;
            if (find_intersection_point.x >= v0.x && find_intersection_point.x <= v1.x &&
                    find_intersection_point.y >= v0.y && find_intersection_point.y <= v3.y)
            {
                tMin = min(tMin, tBase);
            }
        }

        // Check find_intersection with the four triangular faces of the pyramid
        VectorPoint points[4] = {v0, v1, v2, v3};
        for (int i = 0; i < 4; i++)
        {
            VectorPoint v1 = points[i];
            VectorPoint v2 = v4;
            VectorPoint v3 = points[(i + 1) % 4];

            VectorPoint edge1 = v2 - v1;
            VectorPoint edge2 = v3 - v1;
            VectorPoint h = ray.dir.cross(edge2);
            double a = edge1 ^ h;

            if (a > -ZERO && a < ZERO)
                continue;

            double f = 1.0 / a;
            VectorPoint s = ray.start - v1;
            double u = f * (s ^ h);

            if (u < 0.0 || u > 1.0)
                continue;

            VectorPoint q = s.cross(edge1);
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

    virtual Ray find_normal(VectorPoint intersecPoint, Ray incident_ray)
    {
        // Calculate the normal vector at the find_intersection point on the pyramid
        VectorPoint normal = unit_vector(cross_product(VectorPoint(width, 0, 0), VectorPoint(0, width, 0)));
        return Ray(intersecPoint, normal);
    }


    virtual double find_intersect(Ray ray, double *col, int level)
    {
        double tMin = find_tMin(ray);
        if (level == 0)
            return tMin;

        // Check if there was an find_intersection with the pyramid
        if (tMin < INFINITY)
        {
            // Calculate the find_intersection point
            VectorPoint find_intersection_point = ray.start + ray.dir * tMin;

            // Calculate the normal vector at the find_intersection point
            VectorPoint normal;

            // Calculate the normal based on the face of the pyramid that was hit
            VectorPoint apex = lowestPoint + VectorPoint(width / 2, width / 2, height);
            VectorPoint base1 = lowestPoint;
            VectorPoint base2 = lowestPoint + VectorPoint(width, 0, 0);

            if (find_intersection_point.z > apex.z - ZERO)
            {
                // The find_intersection point is on the base of the pyramid
                normal = unit_vector(cross_product(base2 - base1, VectorPoint(0, 0, -1)));
            }
            else
            {
                // The find_intersection point is on one of the triangular faces
                VectorPoint base_mid = lowestPoint + VectorPoint(width / 2, width / 2, 0);

                if (find_intersection_point == apex)
                {
                    // Special case: find_intersection point is at the apex
                    normal = unit_vector(apex - base_mid);
                }
                else
                {
                    // Calculate the normal by taking the cross product of two edges of the triangular face
                    VectorPoint edge1, edge2;
                    if (find_intersection_point == base1)
                    {
                        edge1 = base2 - apex;
                        edge2 = find_intersection_point - apex;
                    }
                    else if (find_intersection_point == base2)
                    {
                        edge1 = base1 - apex;
                        edge2 = find_intersection_point - apex;
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
            VectorPoint reflected_dir = ray.dir - (ray.dir * normal) * normal * 2.0;

            // Set the color of the pyramid at the find_intersection point
            col[0] = color[0] * coEfficients[0];
            col[1] = color[1] * coEfficients[0];
            col[2] = color[2] * coEfficients[0];

            // Recursive reflection
            Ray reflected_ray(find_intersection_point, reflected_dir);
            double reflected_color[3] = {0.0, 0.0, 0.0};
            double tReflect = find_intersect(reflected_ray, reflected_color, level - 1);

            // Add the reflected color to the current color
            col[0] += 0.2 * reflected_color[0];
            col[1] += 0.2 * reflected_color[1];
            col[2] += 0.2 * reflected_color[2];

            return tMin;
        }
        return -1.0; // No find_intersection
    }

    void setLowestPoint(const VectorPoint& point)
    {
        lowestPoint = point;
    }

    VectorPoint getLowestPoint() const
    {
        return lowestPoint;
    }

    void setWidth(double w)
    {
        width = w;
    }

    double getWidth() const
    {
        return width;
    }

    void setHeight(double h)
    {
        height = h;
    }

    double getHeight() const
    {
        return height;
    }

    void setReflectivity(double reflectivity)
    {
        // Set reflectivity
    }

    double getReflectivity() const
    {
        // Get reflectivity
    }

    void setTransparency(double transparency)
    {
        // Set transparency
    }

    double getTransparency() const
    {
        // Get transparency
    }

};

void Pyramid::input_pyramid(ifstream &ifs)
{
    ifs >> lowestPoint >> width >> height;
    input_obj(ifs);
}

void Pyramid::print()
{
    print_object_details();
}


