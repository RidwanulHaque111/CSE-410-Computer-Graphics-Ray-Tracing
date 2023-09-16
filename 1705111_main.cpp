#include <GL/glut.h>

#include <bits/stdc++.h>
#include "1705111_objects.hpp"
#include "bitmap_image.hpp"


//Global data
int imageWidth, imageHeight;
int consoleWidth = 600, consoleHeight = 600;
int floor_shine = 0;
int image_count = 0;


double nearPlaneDistance, farPlaneDistance, fovY, aspectRatio;
double floor_amb, floor_diff, floor_refl;
double floor_spec = 0.0;
int tileSize;

int drawgrid;
int drawaxes;


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void photoCapture()
{
    cout << "\n----------RAY TRACING LOG----------\n";

	// Initialize image
	bitmap_image image(imageWidth, imageHeight);
	for (int i = 0; i < imageWidth; i++)
		for (int j = 0; j < imageHeight; j++)
			image.set_pixel(i, j, 0, 0, 0);

    double planeDistance = (consoleHeight/2.0) / tan(rad(fovY/2.0));
    VectorPoint topLeft = camera_pos + l*planeDistance - r*(consoleWidth/2.0) + u*(consoleHeight/2.0);

    double du = (1.0*consoleWidth)/(1.0*imageWidth);
    double dv = (1.0*consoleHeight)/(1.0*imageHeight);

    double *color = new double[3];
    VectorPoint tl_pixel = topLeft + r*(du*0.5) - u*(dv*0.5);

    int nearest;
    double t, tMin;
    for(int i=0; i<imageWidth; i++) {
        for(int j=0; j<imageHeight; j++) {
            VectorPoint current_pixel = tl_pixel + r*(i*du) - u*(j*dv);
            Ray ray(camera_pos, current_pixel-camera_pos);
            VectorPoint startPoint = ray.getStartPoint();
            VectorPoint endPoint = ray.getDirPoint();

            //Closest object
            nearest = INT_MAX;
            tMin=INT_MAX;
            for(int k=0; k<objects.size(); k++) {
                color[0] = 0.0; color[1] = 0.0; color[2] = 0.0;
                t = objects.at(k)->find_intersect(ray, color, 0);
                if(t>0.0 && t<tMin) {
                    tMin = t; nearest = k;
                }
            }

            //Assign color of closest object (if it exists)
            if(nearest < INT_MAX) {
                color[0] = 0.0; color[1] = 0.0; color[2] = 0.0;
                tMin = objects[nearest]->find_intersect(ray, color, 1);
                color[0] *= 255.0; color[1] *= 255.0; color[2] *= 255.0;
                image.set_pixel(i, j, round(color[0]), round(color[1]), round(color[2]));
            }
        }
    }

    delete[] color;
	cout << "Preparing the Image ......\n\n";
    string image_name = "D:\\BUET_Files\\L4T1\\CSE_410_Computer_Graphics\\offline3\\ray_tracing\\output\\out_" + std::to_string(++image_count) + ".bmp";
	image.save_image(image_name);
    cout <<"Success : "<< "Image "<<"\"out_"<<image_count<<"\" Saved.\n\n";
}


template<typename Iterator>
void deleteObjects(Iterator start, Iterator last)
{
   while ( start != last ){
       delete *start;
       start++;
   }
}

//dealloc the vectors
void clearMemory()
{
    cout << "Clearing Memory......\n\n";
    deleteObjects(objects.begin(), objects.end());
    objects.clear();
    cout << "All Objects Cleared ...\n";
    deleteObjects(lights_list.begin(), lights_list.end());
    lights_list.clear();
    cout << "All Lights Cleared ...\n\n";
    cout << "Successfully Cleared All Memory.\n";
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){
				drawaxes=1-drawaxes;
			}
			break;
		case GLUT_RIGHT_BUTTON:
			//........
			break;
		case GLUT_MIDDLE_BUTTON:
			//........
			break;
		default:
			break;
	}
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			//drawgrid=1-drawgrid;
			axis_rotation(&r, &l, &u, 3);
			break;

		case '2':
			axis_rotation(&r, &l, &u, -3);
			break;

		case '3':
			axis_rotation(&l, &u, &r, 3);
			break;

		case '4':
			axis_rotation(&l, &u, &r, -3);
			break;

		case '5':
			axis_rotation(&u, &r, &l, 3);
			break;

		case '6':
			axis_rotation(&u, &r, &l, -3);
			break;

		case '0':
			photoCapture();
			break;

		case 'q':
            exit(0);
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			camera_pos = move_along_unit_vect(camera_pos, l, -2);
			break;
		case GLUT_KEY_UP:		// up arrow key
			camera_pos = move_along_unit_vect(camera_pos, l, 2);
			break;

		case GLUT_KEY_RIGHT:
			camera_pos = move_along_unit_vect(camera_pos, r, 2);
			break;
		case GLUT_KEY_LEFT:
			camera_pos = move_along_unit_vect(camera_pos, r, -2);
			break;

		case GLUT_KEY_PAGE_UP:
            camera_pos = move_along_unit_vect(camera_pos, u, 2);
			break;
		case GLUT_KEY_PAGE_DOWN:
            camera_pos = move_along_unit_vect(camera_pos, u, -2);
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;

		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    gluLookAt(camera_pos.x,camera_pos.y,camera_pos.z,    camera_pos.x+l.x,camera_pos.y+l.y,camera_pos.z+l.z,   u.x,u.y,u.z);

	glMatrixMode(GL_MODELVIEW);

	//drawing in OpenGL
	drawAxes();
    for(int i=0; i<objects.size(); i++)
        objects.at(i)->draw();

    for(int i=0; i<lights_list.size(); i++)
        lights_list.at(i)->draw();

	glutSwapBuffers();
}


void animate(){
	glutPostRedisplay();
}

void init(){
	drawgrid=0;
	drawaxes=0;

	//camera vectors
    u.x = 0; u.y = 0; u.z = 1;
    r.x = -1/sqrt(2); r.y = 1/sqrt(2); r.z = 0;
    l.x = -1/sqrt(2); l.y = -1/sqrt(2); l.z = 0;
    camera_pos.x = 80.0; camera_pos.y = 80.5; camera_pos.z = 22.5;

	//clear the screen
	glClearColor(0,0,0,0);

	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fovY, aspectRatio, nearPlaneDistance, farPlaneDistance);
}


void runRTX(){

    ifstream inFile("D:\\BUET_Files\\L4T1\\CSE_410_Computer_Graphics\\offline3\\ray_tracing\\scene.txt");
    if(inFile.is_open()){
        inFile >> nearPlaneDistance >> farPlaneDistance >> fovY >> aspectRatio;
        inFile >> recursion_level >> imageHeight;
        imageWidth = imageHeight * aspectRatio;

        inFile >> tileSize;
        Floor *f = new Floor(farPlaneDistance, tileSize);
        inFile >> floor_amb >> floor_diff >> floor_refl;
        f->setColors(1.0, 1.0, 1.0);
        f->setCoEfficients(floor_amb, floor_diff, floor_spec, floor_refl);
        f->setShine(floor_shine);
        objects.push_back(f);

        int n; string name;
        inFile >> n;

        //reading objects
        for(int i=0; i<n; i++){
            inFile >> name;
            if(!name.compare("sphere")) {
                Sphere *sp = new Sphere();
                sp->input_sphere(inFile);
                objects.push_back(sp);
            }
			else if(!name.compare("pyramid")) {
				Pyramid *py = new Pyramid();
				py->input_pyramid(inFile);
				objects.push_back(py);
			}
			else if(!name.compare("cube")) {
				Cube *cb = new Cube();
				cb->input_cube(inFile);
				objects.push_back(cb);
			}
        }

        inFile >> n;
        //scan point lights from file
        for(int i=0; i<n; i++){
            PointLightObj *pl = new PointLightObj();
            pl->input_pointlight(inFile);
            lights_list.push_back(pl);
        }

        inFile >> n;
        //scan spot lights from file
        for(int i=0; i<n; i++){
            SpotLight *sl = new SpotLight();
            sl->input_spotlight(inFile);
            lights_list.push_back(sl);
            cout<<"\nAll ready to Go!"<<endl;
            cout<<"----------------"<<endl;
        }
    } else {
        cout << "Please Enter a correct File!" << endl;
        exit(0);
    }
}


int main(int argc, char **argv){
    //input from scene.txt
    std::cout << std::setprecision(7) << std::fixed;
    atexit(clearMemory);
    runRTX();

	glutInit(&argc,argv);
	glutInitWindowSize(consoleHeight, consoleWidth);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");
	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
