class Sphere {
	private:
	/***Characteristics***/
	double radius;
	double x_origin;
	double y_origin;
	double z_origin;
	double originNorm;
	
	public:
	/***Constructors***/
	Sphere(double r, double x, double y, double z) {
		radius = r;
		x_origin = x;
		y_origin = y;
		z_origin = z;
	}
	Sphere(double r) {
		radius = r;
		x_origin = 0;
		y_origin = 0;
		z_origin = 0;
	}
	
	void setOriginNorm(double o) { originNorm = o; }
	
	/***Methods***/
	double getRadius() { return radius; }
	double getOrigin_X() { return x_origin; }
	double getOrigin_Y() { return y_origin; }
	double getOrigin_Z() { return z_origin; }
	double getOriginNorm() { return originNorm; }
};