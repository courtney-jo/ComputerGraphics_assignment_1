// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include "sphere.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;

/*
	//calculate shade ray
	Vector3d s = ray_intersection + (t*l);
	
	//ray blocked by another object
	if(s.norm() >= 0){
		Vector3d h = l + (ray_direction*-1);
		return h.normalized().transpose() * ray_normal; 
	}*/

/***Extra Functions***/
double diffuse(const Vector3d lightPos, Vector3d ray_intersection) {
	//compute intersection point and normal
	Vector3d l = lightPos - ray_intersection;
    Vector3d ray_normal = ray_intersection.normalized();

	return l.normalized().transpose() * ray_normal;
}

double specular(const Vector3d lightPos, Vector3d ray_intersection,  Vector3d origin) {
	//Compute intersection point and normal
    Vector3d ray_normal = ray_intersection.normalized();

	//Specular model
	Vector3d h = (origin-ray_intersection).normalized().transpose() + (lightPos - ray_intersection).normalized().transpose();
	return pow(ray_normal.dot(h/h.norm()),100);
}

double hit(Vector3d s_origin, Vector3d ray_origin, Vector3d ray_direction, double r) {
	Vector3d p = ray_origin - s_origin;
	double t1 = -1 * ray_direction.dot(p);
	double t2 = ray_direction.dot(ray_direction) * (p.dot(p) - (r*r));
	double t3 = pow(ray_direction.dot(p),2);
	double t4 = sqrt(t3 - t2);
	double t = (t1+t4) / ray_direction.dot(ray_direction);
	double other_t = (t1-t4) / ray_direction.dot(ray_direction);
	
	if(t < other_t)
		return t;
	else
		return other_t;
}

/***Assignment Parts***/

void part1_1()
{
    std::cout << "Part 1_1: multiple spheres with different origins" << std::endl;

    const std::string filename("part1_1.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/C.cols(),0,0);
    Vector3d y_displacement(0,-2.0/C.rows(),0);

    // Single light source
    const Vector3d light_position(-1,1,1);
	
	/***Spheres***/
	Sphere S1(0.5,.2,.1,.1);
	Vector3d s1_origin(S1.getOrigin_X(),S1.getOrigin_Y(),S1.getOrigin_Z());
	
	Sphere S2(0.2,-.2,0,.4);
	Vector3d s2_origin(S2.getOrigin_X(),S2.getOrigin_Y(),S2.getOrigin_Z());

	
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

            // Intersect with the sphere
			double t1 = hit(s1_origin,ray_origin,ray_direction, S1.getRadius());
			double t2 = hit(s2_origin,ray_origin,ray_direction, S2.getRadius());
			
			/***S1***/
            if (t1 > 0)
            {
				Vector3d ray_intersection=ray_origin+(t1*ray_direction);
				
				C(i,j) = diffuse(light_position, ray_intersection);
                C(i,j) = max(C(i,j),0.);	//color and clamp to 0
                A(i,j) = 1;	//disable alpha mask
            }
			
            /***S2***/
            if (t2 > 0)
            {
				Vector3d ray_intersection=ray_origin+(t2*ray_direction);
				
				C(i,j) = diffuse(light_position, ray_intersection);
                C(i,j) = max(C(i,j),0.);	//color and clamp to 0
                A(i,j) = 1;	//disable alpha mask
            }
        }
    }

    // Save to png
    write_matrix_to_png(C,C,C,A,filename);
}

void part1_2()
{
    std::cout << "Part 1_2:  multiple spheres, different colors, different materials, two lights" << std::endl;

    const std::string filename("part1_2.png");
	//Size of Matrix determines camera view size
    MatrixXd R = MatrixXd::Zero(800,800);	//Store Red intensity
	MatrixXd G = MatrixXd::Zero(800,800);	//Store Green intensity
	MatrixXd B = MatrixXd::Zero(800,800);	//Store Blue intensity
    MatrixXd A = MatrixXd::Zero(800,800);	// Store the alpha mask
	
	//Materials
	MatrixXd D = MatrixXd::Zero(800,800);	MatrixXd D2 = MatrixXd::Zero(800,800);
	MatrixXd S = MatrixXd::Zero(800,800);	MatrixXd S2 = MatrixXd::Zero(800,800);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    //Set light source
    const Vector3d light_position(-1,1,1);
	const Vector3d light_position2(1,1,1);
	
	/***Spheres***/
	Sphere S1(0.4,-.3,-.2,-.3);
	Vector3d s1_origin(S1.getOrigin_X(),S1.getOrigin_Y(),S1.getOrigin_Z());
	
	Sphere Sp2(0.4,.2,.1,.2);
	Vector3d s2_origin(Sp2.getOrigin_X(),Sp2.getOrigin_Y(),Sp2.getOrigin_Z());

    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {		
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

			//ray intersection
			double t1 = hit(s1_origin,ray_origin,ray_direction, S1.getRadius());
			double t2 = hit(s2_origin,ray_origin,ray_direction, Sp2.getRadius());
			
			/***S1***/
            if (t1 > 0)
            {
				Vector3d ray_intersection=ray_origin+(t1*ray_direction);
				
				D(i,j) = diffuse(light_position, ray_intersection);
				D2(i,j) = diffuse(light_position2, ray_intersection);

                // Clamp to zero and set color
                R(i,j) = .5 * (.3 + max(D(i,j),0.) + max(D2(i,j),0.));
				G(i,j) = 0 * (max(D(i,j),0.) + max(D2(i,j),0.));
				B(i,j) = .5 * (.3 + max(D(i,j),0.) + max(D2(i,j),0.));

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
			
            /***Sp2***/
            if (t2 > 0)
            {
				Vector3d ray_intersection=ray_origin+(t2*ray_direction);
				
				// diffuse
				D(i,j) = diffuse(light_position, ray_intersection);
				D2(i,j) = diffuse(light_position2, ray_intersection);
			
                // specular
				S(i,j) = specular(light_position, ray_intersection, origin);
				S2(i,j) = specular(light_position2, ray_intersection, origin);

                // Clamp to zero and set color
                R(i,j) = 0 * (max(D(i,j),0.) + max(S(i,j),0.) + max(D2(i,j),0.) + max(S2(i,j),0.));
				G(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.) + max(D2(i,j),0.) + max(S2(i,j),0.));
				B(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.) + max(D2(i,j),0.) + max(S2(i,j),0.));

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);
}

void part1_3_1()
{
    std::cout << "Part 1_3_1: perspective projection" << std::endl;

    const std::string filename("part1_3_1.png");
	//Size of Matrix determines camera view size
    MatrixXd C = MatrixXd::Zero(800,800);	//Store color
    MatrixXd A = MatrixXd::Zero(800,800);	// Store the alpha mask

    //camera size and angle
    Vector3d origin(1,-1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    //Set light source
    const Vector3d light_position(1,-1.15,.2);
	
	/***Spheres***/
	Sphere S1(0.5,1.2,-1,-.4);
	Vector3d s1_origin(S1.getOrigin_X(),S1.getOrigin_Y(),S1.getOrigin_Z());
	
	Sphere S2(0.2,.8,-1.15,0);
	Vector3d s2_origin(S2.getOrigin_X(),S2.getOrigin_Y(),S2.getOrigin_Z());

    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {	
			/***Perspective Projection***/
			
			//Prepare the ray
			Vector3d ray_direction =  double(i)*x_displacement + double(j)*y_displacement - origin;
			ray_direction.normalize();
			
			// Intersect with the sphere
			double t1 = hit(s1_origin,origin,ray_direction, S1.getRadius());
			double t2 = hit(s2_origin,origin,ray_direction, S2.getRadius());
			
			/***S1***/
            if (t1 > 0)
            {
				Vector3d ray_intersection=origin+(t1*ray_direction);
				
				C(i,j) = diffuse(light_position, ray_intersection);
                C(i,j) = max(C(i,j),0.);	//color and clamp to zero
                A(i,j) = 1;	// Disable the alpha mask for this pixel
            }
			
            /***S2***/
            if (t2 > 0)
            {
				Vector3d ray_intersection=origin+(t2*ray_direction);
				
				C(i,j) = diffuse(light_position, ray_intersection);
                C(i,j) = max(C(i,j),0.);	//color and clamp to zero
                A(i,j) = 1;	//disable alpha mask
            }
        }
    }

    // Save to png
    write_matrix_to_png(C,C,C,A,filename);
}

void part1_3_2()
{
	/*perspective good make tiny adjustments to placement and fix lighting*/
    std::cout << "Part 1_3_2: perspective projection" << std::endl;

    const std::string filename("part1_3_2.png");
	//Size of Matrix determines camera view size
    MatrixXd R = MatrixXd::Zero(800,800);	//Store Red intensity
	MatrixXd G = MatrixXd::Zero(800,800);	//Store Green intensity
	MatrixXd B = MatrixXd::Zero(800,800);	//Store Blue intensity
    MatrixXd A = MatrixXd::Zero(800,800);	// Store the alpha mask
	
	//Materials
	MatrixXd D = MatrixXd::Zero(800,800);	MatrixXd D2 = MatrixXd::Zero(800,800);
	MatrixXd S = MatrixXd::Zero(800,800);	MatrixXd S2 = MatrixXd::Zero(800,800);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(1,-1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    //Set light source
    const Vector3d light_position(1,-1,.5);
	const Vector3d light_position2(-1,1,.5);
	
	/***Spheres***/
	Sphere S1(0.4,1.2,-1,-.2);
	Vector3d s1_origin(S1.getOrigin_X(),S1.getOrigin_Y(),S1.getOrigin_Z());
	
	Sphere Sp2(0.4,.6,-.7,.0);
	Vector3d s2_origin(Sp2.getOrigin_X(),Sp2.getOrigin_Y(),Sp2.getOrigin_Z());

    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {		
			/***Perspective Projection***/
			
			//Prepare the ray
			Vector3d ray_direction =  double(i)*x_displacement + double(j)*y_displacement - origin;
			ray_direction.normalize();
			
			
			// Intersect with the sphere
			double t1 = hit(s1_origin,origin,ray_direction, S1.getRadius());
			double t2 = hit(s2_origin,origin,ray_direction, Sp2.getRadius());
			
			/***S1***/
            if (t1 > 0)
            {
				Vector3d ray_intersection=origin+(t1*ray_direction);
				
				D(i,j) = diffuse(light_position, ray_intersection);
				D2(i,j) = diffuse(light_position2, ray_intersection);

                // Clamp to zero and set color
                R(i,j) = .5 * (.3 + max(D(i,j),0.) + max(D2(i,j),0.));
				G(i,j) = 0 * (max(D(i,j),0.) + max(D2(i,j),0.));
				B(i,j) = .5 * (.3 + max(D(i,j),0.) + max(D2(i,j),0.));

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
			
            /***Sp2***/
            if (t2 > 0)
            {
				Vector3d ray_intersection=origin+(t2*ray_direction);
				
				// diffuse
				D(i,j) = diffuse(light_position, ray_intersection);
				D2(i,j) = diffuse(light_position2, ray_intersection);
			
                // specular
				S(i,j) = specular(light_position, ray_intersection, origin);
				S2(i,j) = specular(light_position2, ray_intersection, origin);

                // Clamp to zero and set color
                R(i,j) = 0 * (max(D(i,j),0.) + max(S(i,j),0.) + max(D2(i,j),0.) + max(S2(i,j),0.));
				G(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.) + max(D2(i,j),0.) + max(S2(i,j),0.));
				B(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.) + max(D2(i,j),0.) + max(S2(i,j),0.));

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);
}

void part1_5()
{
    std::cout << "Part 1_5: shadows" << std::endl;

    const std::string filename("part1_5.png");
	//Size of Matrix determines camera view size
    MatrixXd R = MatrixXd::Zero(800,800);	//Store Red intensity
	MatrixXd G = MatrixXd::Zero(800,800);	//Store Green intensity
	MatrixXd B = MatrixXd::Zero(800,800);	//Store Blue intensity
    MatrixXd A = MatrixXd::Zero(800,800);	// Store the alpha mask
	
	//Materials
	MatrixXd D = MatrixXd::Zero(800,800);
	MatrixXd S = MatrixXd::Zero(800,800);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    //Set light source
    const Vector3d light_position(-.6,.6,.4);
	
	/***Spheres***/
	Sphere S1(0.5,.6,.2,0);
	Vector3d s1_origin(S1.getOrigin_X(),S1.getOrigin_Y(),S1.getOrigin_Z());
	
	Sphere S2(0.35,-.3,.2,0);
	Vector3d s2_origin(S2.getOrigin_X(),S2.getOrigin_Y(),S2.getOrigin_Z());

    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {		
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

			//ray intersection
			double t1 = hit(s1_origin,ray_origin,ray_direction, S1.getRadius());
			double t2 = hit(s2_origin,ray_origin,ray_direction, S2.getRadius());
			
			//object hit
			if(t1 > 0 || t2 > 0) {
				/***S1***/
				if (t1 > 0)
				{
					Vector3d ray_intersection=ray_origin+(t1*ray_direction);
						
					// diffuse
					D(i,j) = diffuse(light_position, ray_intersection);
					
					// specular
					S(i,j) = specular(light_position, ray_intersection, origin);

					// Clamp to zero and set color
					R(i,j) = 0 * (max(D(i,j),0.) + max(S(i,j),0.));
					G(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.));
					B(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.));

					// Disable the alpha mask for this pixel
					A(i,j) = 1;
					
					//shadow
					//create ray going towards other object from the point we just hit and see if we hit the object
					//looks like it works but doesn't
					double shadow_ray = hit(s2_origin,ray_origin-s1_origin,ray_direction,S1.getRadius());
					if (shadow_ray >= 0) {
						R(i,j) = 0 * 0;
						G(i,j) = .5 * .3;
						B(i,j) = .5 * .3;

						// Disable the alpha mask for this pixel
						A(i,j) = 1;
					}
				}

				/***S2***/
				if (t2 > 0)
				{
					Vector3d ray_intersection=ray_origin+(t2*ray_direction);
					
					D(i,j) = diffuse(light_position, ray_intersection);

					// Clamp to zero and set color
					R(i,j) = .5 * (.3 + max(D(i,j),0.));
					G(i,j) = 0 * (max(D(i,j),0.));
					B(i,j) = .5 * (.3 + max(D(i,j),0.));

					// Disable the alpha mask for this pixel
					A(i,j) = 1;
					
					double shadow_ray = hit(s1_origin,ray_origin-s2_origin,ray_direction,S2.getRadius());
					if (shadow_ray >= 0) {
						R(i,j) = .5 * .3;
						G(i,j) = 0 * 0;
						B(i,j) = .5 * .3;

						// Disable the alpha mask for this pixel
						A(i,j) = 1;
					}
				}
			}
			//background
			else {
				// Clamp to zero and set color
				R(i,j) = .5;
				G(i,j) = .5;
				B(i,j) = .5;

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
				
				double shadow_ray = hit(ray_origin,ray_origin-s2_origin,ray_direction,S2.getRadius());
					if (shadow_ray >= 0) {
						R(i,j) = 0;
						G(i,j) = 0;
						B(i,j) = 0;

						// Disable the alpha mask for this pixel
						A(i,j) = 1;
					}
					shadow_ray = hit(ray_origin,ray_origin-s1_origin,ray_direction,S1.getRadius());
					if (shadow_ray >= 0) {
						R(i,j) = 0;
						G(i,j) = 0;
						B(i,j) = 0;

						// Disable the alpha mask for this pixel
						A(i,j) = 1;
					}
			}
        }
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);
}

void part1_6()
{
    std::cout << "Part 1_6: reflection" << std::endl;

    const std::string filename("part1_6.png");
	//Size of Matrix determines camera view size
    MatrixXd R = MatrixXd::Zero(800,800);	//Store Red intensity
	MatrixXd G = MatrixXd::Zero(800,800);	//Store Green intensity
	MatrixXd B = MatrixXd::Zero(800,800);	//Store Blue intensity
    MatrixXd A = MatrixXd::Zero(800,800);	// Store the alpha mask
	
	//Materials
	MatrixXd D = MatrixXd::Zero(800,800);
	MatrixXd S = MatrixXd::Zero(800,800);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    //Set light source
    const Vector3d light_position(-1,1,1);
	
	/***Spheres***/
	Sphere S1(0.35,-.32,.2,0);
	Vector3d s1_origin(S1.getOrigin_X(),S1.getOrigin_Y(),S1.getOrigin_Z());
	
	Sphere S2(0.35,.32,.2,0);
	Vector3d s2_origin(S2.getOrigin_X(),S2.getOrigin_Y(),S2.getOrigin_Z());

    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {		
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

			//ray intersection
			double t1 = hit(s1_origin,ray_origin,ray_direction, S1.getRadius());
			double t2 = hit(s2_origin,ray_origin,ray_direction, S2.getRadius());
			
			//object hit
			if(t1 > 0 || t2 > 0) {
				/***S1***/
				if (t1 > 0)
				{
					Vector3d ray_intersection=ray_origin+(t1*ray_direction);
						
					// diffuse
					D(i,j) = diffuse(light_position, ray_intersection);
					
					// specular
					S(i,j) = specular(light_position, ray_intersection, origin);

					// Clamp to zero and set color
					R(i,j) = 0 * (max(D(i,j),0.) + max(S(i,j),0.));
					G(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.));
					B(i,j) = .5 * (.3 + max(D(i,j),0.) + max(S(i,j),0.));
					
					// Disable the alpha mask for this pixel
					A(i,j) = 1;
					
					double shadow_ray = hit(s2_origin,ray_origin-s1_origin,ray_direction,S1.getRadius());
					if (shadow_ray >= 0) {
						R(i,j) = 0;
						G(i,j) = 0;
						B(i,j) = 0;

						// Disable the alpha mask for this pixel
						A(i,j) = 1;
					}
				}
				
				/***S2***/
				if (t2 > 0)
				{
					Vector3d ray_intersection=ray_origin+(t2*ray_direction);
					
					D(i,j) = diffuse(light_position, ray_intersection);

					// Clamp to zero and set color
					R(i,j) = .5 * (.3 + max(D(i,j),0.));
					G(i,j) = 0 * (max(D(i,j),0.));
					B(i,j) = .5 * (.3 + max(D(i,j),0.));

					// Disable the alpha mask for this pixel
					A(i,j) = 1;
				}
			}
			//background
			else if (ray_origin(1) < .3) {
				//floor for reflection
				//would add the recursive color function here along with mirror material
				// Clamp to zero and set color
				R(i,j) = 1;
				G(i,j) = 1;
				B(i,j) = 1;

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
			else {
				// Clamp to zero and set color
				R(i,j) = .5;
				G(i,j) = .5;
				B(i,j) = .5;

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
        }
    }

    // Save to png
    write_matrix_to_png(R,G,B,A,filename);
}

int main()
{
	part1_1();		//complete have two spheres of different origins
	part1_2();		//two spheres, different colors, different materials
	part1_3_1();	//prespective projection
	part1_3_2();	//prespective projection
	part1_5();		//shadows
	part1_6();	//reflection

    return 0;
}
