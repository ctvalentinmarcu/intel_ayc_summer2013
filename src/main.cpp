/*!
 * \file main.cpp
 * \author Cédric Andreolli (Intel)
 * \date 10 April 2013
 * This file contains the main algorithm used to match templates in a bmp file. The algorithm here is very simple and completly not optimized. 
 * Here are different things you can do:
 *<ul><li>Improve the algorithm:
 *		<ul><li>Take the rotations in consideration</li>
 * 		<li>Take a finesse granularity of scales</li>
 * 		<li>Take lights and luminosity in consideration</li>
 *		<li>Your goal is to create a good and fast algorithm to retrieve a pattern in a BMP file</li></ul>
 *   </li>
 *   <li>Optimize your algorithm:
 * 		<ul><li>Be careful with the memory utilization</li>
 *		<li>Avoid cache misses</li>
 *		<li>Use vectorization and Intrinsic functions</li></ul>
 *  </li>
 *<ul>
 * The notation will take in consideration your results on the xeon phi. Programming for the xeon phi is not easy (even if a classic C++ code can run natively on this card),
 * you will need to think your program as a program able to scale on massive parallel hardware. To help you debug and optimize your program, feel free to use Intel softwares
 * such as VTune Amplifier or Inspector (free one year licenses are available on <a href="http://www.intel-software-academic-program.com">this website</a>.
 */


/*! \mainpage Accelerate Your Code contest - Early 2013 
 *
 * \section intro Introduction
 * This section describes the code sample provided for the Intel* Accelerate Contests 2013 - early edition.
 * \section todo What do you have to do?
 * The code provided is just a very simple implementation of a pattern matching algorithm. You need to modify and improve this algorithm to achieve 2 goals:
 * <ul>
 *	<li>Create a better algorithm</li>
 *	<li>Create a faster algorithm</li>
 * </ul>
 * Those 2 aspects will be evaluated by Intel engineers and will be taken in consideration to choose the winners.
 * \section better_algo Create a better algorithm
 * By creating a better algorithm, we mean that your algorithm should be able to introduce rotation, luminosity and what ever you think important in a good pattern matching algorithm. You can also play with the scale factor
 * with more finesse. 
 * \section faster_algo Create a faster algorithm
 * To make your algorithm faster, you will need to optimize it. Optimization can take a lot of different aspects:
 * <ul>
 * 	<li>Avoiding cache misses</li>
 *	<li>Using compiler options</li>
 *	<li>Using vectorization</li>
 *	<li>Working on the algorithm to avoid un-necessarry computations</li>
 * </ul>
 * We will provide a simple documentation about this specific topic during the contest. 
 * \section xeon_phi The Xeon Phi
 * For this contest, you will have the chance to use the new Intel Xeon Phi. The Xeon Phi is a co-processor plugged on the motherboard via PCI. It allows to run 240 threads at a time. The great thing with Xeon Phi is that it is
 * able to run C++ native code. You also have access to 512 bits registers that gives vectorization a new meaning :)
 * But programming for such a device means that your program must be very well optimized to handle such a massive parallelism. A special attention will be given to the candidates who have good performances on this device.
 * \section questions Questions?
 * If you have any questions feel free to contact us though the forum or directly by email:
 * <ul>
 *	<li><a href="mailto:paul.guermonprez@intel.com" >Paul Guermonprez</a> - Academic Program Manager</li>
 *	<li><a href="mailto:cedric.andreolli@intel.com" >Cédric Andreolli</a> - Academic Program Software Developer Intern</li>
 * </ul>
 */

#include <iostream>
#include "BitmapImporter.h"
#include <omp.h>
#include <string>
#include <list>
#include<vector>
#include <stdlib.h>
#include <iterator>

using namespace std;
using namespace Accelerate;

unsigned char maxdif = 10;
float maxwrong_ratio = 0.3;
float maxwrong_ratio_edge = 0.6;

/*!
 *\struct Parameters
 *\brief This structure holds the application parameters
 */
typedef struct{
	int nb_threads;
	string main_image_name;
	list<string> template_names;
	int max_scale;
}Parameters;

/*!
 *\brief This structure hold a single result
 */
typedef struct{
	int pattern_ID;
	int position_x;
	int position_y;
}Result;


void calculate_transform(Affine* transform, const Point& p1, const Point& p2, const Point& p3, const Point& p4);

bool match_point(const Point& p1, const Point& p2);

Point apply_transform(const Point& p, const Affine& transform);


/*!
 * \brief Try to match the exact template in the main image starting at coordinates [h,w] in the main image.
 * \param main_image The main image.
 * \param template_image The template.
 * \param h The row number.
 * \param w the coloumn number.
 * \return template has been found at position [h,w] in the main image.  
 */
bool match_template(const Accelerate::Image& main_image, const Accelerate::Image& template_image, unsigned int &h, unsigned int &w, int max_scale);

bool match1(const Accelerate::Image& main_image, const Accelerate::Image& template_image, unsigned int &h, unsigned int &w, int max_scale);

/*!
 * \brief Read the parameters
 * \param argc The number of parameters (including the program call).
 * \param argv The array of parameters.
 * \return The parameters have been read correctly.
 */
bool read_parameters(int argc, char* argv[], Parameters& parameters);


bool compare_results(Result& first, Result& second){
	if(first.pattern_ID < second.pattern_ID) return true;
	else if(first.pattern_ID > second.pattern_ID) return false;

	if(first.position_x < second.position_x) return true;
	else if(first.position_x > second.position_x) return false;

	if(first.position_y < second.position_y) return true;
	else if(first.position_y > second.position_y) return false;

	return true;	
}

bool read_parameters(int argc, char* argv[], Parameters& parameters){
	if(argc < 4) return false;

	parameters.nb_threads = atoi(argv[1]);
	if(parameters.nb_threads < 0) return false;

	parameters.max_scale = atoi(argv[2]);
	if(parameters.max_scale <= 0) return false;

	parameters.main_image_name = string(argv[3]);

	for(unsigned int i=4; i<argc; i++){
		parameters.template_names.push_back(string(argv[i]));
	}
	return true;	
}

int main(int argc, char* argv[]){
	if(argc<=1) return 0;
	if(argv == NULL) return 0;

	Parameters parameters;
	list<Result> result_list;
	if(!read_parameters(argc, argv, parameters)){
		cout<<"Wrong number of parameters or invalid parameters..."<<endl;
		cout<<"The program must be called with the following parameters:"<<endl;
		cout<<"\t- num_threads: The number of threads"<<endl;
		cout<<"\t- max_scale: The maximum scale that can be applied to the templates in the main image"<<endl;
		cout<<"\t- main_image: The main image path"<<endl;
		cout<<"\t- t1 ... tn: The list of the template paths. Each template separated by a space"<<endl;
		cout<<endl<<"For example : ./run 4 3 img.bmp template1.bmp template2.bmp"<<endl;
		return -1;
	}
	if(parameters.nb_threads > 0)
		omp_set_num_threads(parameters.nb_threads);
	else
		omp_set_num_threads(240);

	//Read the main image
	Accelerate::Image main_image = Accelerate::Image::create_image_from_bitmap(parameters.main_image_name, parameters.max_scale, 0, 0);

//#pragma omp parallel //for schedule(dynamic)

//#pragma omp single

	//iterates over the pattern images
	for(string temp_name : parameters.template_names){ 
		//Read a specific pattern image
		Accelerate::Image template_image = Accelerate::Image::create_image_from_bitmap(temp_name, parameters.max_scale, main_image.get_height(), main_image.get_width());
		//Then retrieve the ID of the image (the 3 first characters
		int template_id = atoi(temp_name.substr(0, 3).c_str());	
		//Iterate over some possible scales (you can add more steps and you can also check rotations)
		//The sample is really simple because you have to create a good algorithm able to match
		//a pattern in an image
		//for(unsigned int s=1; s<=parameters.max_scale; s++){
			//Create a scaled image
			//Accelerate::Image temp = template_image.scale_image(s);
			//iterates on the main image pixels		

		//Point first_edge = t.reference;
#pragma omp parallel //for schedule(dynamic)
#pragma omp single
		for(unsigned int hm=0; hm<main_image.get_height(); hm++){
			for(unsigned int wm=0; wm<main_image.get_width(); wm++){			
				//Try to match the template
				//#pragma omp task untied firstprivate(wm,hm)
//#if 1
				Point first_edge = template_image.reference;
				PixelStr pixel = main_image.pixel_data[main_image.get_width()*hm + wm];
				PixelStr first = template_image.pixel_data[0];
				if(abs(pixel.r - first_edge.r) > maxdif || abs(pixel.g - first_edge.g) > maxdif || abs(pixel.b - first_edge.b) > maxdif)
					continue;
#if 0				
				bool found = false;						
				for(Point c : main_image.edgeR[template_image.index_min])
					if(match_point(first_edge, c))
						found = true;		
				
				for(Point c : main_image.edgeG[template_image.index_min])
					if(match_point(first_edge, c))
						found = true;					
				for(Point c : main_image.edgeB[template_image.index_min])
					if(match_point(first_edge, c))
						found = true;				
				if(found == false) 
					continue;
#endif
				
				unsigned int reth = hm;
				unsigned int retw = wm;

				#pragma omp task untied firstprivate(retw,reth)
				if(match1(main_image, template_image, reth, retw, parameters.max_scale)){
				//if(match_template(main_image, template_image, reth, retw, parameters.max_scale)){
					//The pattern image has been found so save the result
					Result result;
					result.pattern_ID = template_id;
					result.position_x = retw;
					result.position_y = reth;
					#pragma omp critical
					result_list.push_back(result);		
				}
			}
		}

	
		//}
	}
	//sort the result list
	result_list.sort(compare_results);
	int lastx = -100;
	int lasty = -100;
	int mindif = 10;
	//Print the results
	for(Result res : result_list) 
		if(abs(res.position_x - lastx) > mindif && abs(res.position_y - lasty) > mindif)
		{
			lastx = res.position_x;
			lasty = res.position_y;
			cout<<res.pattern_ID<<"\t"<<res.position_x<<"\t"<<res.position_y<<endl;
		}

	return 0;
}


bool match1(const Accelerate::Image& main, const Accelerate::Image& t, unsigned int &h, unsigned int &w, int max_scale){
	//unsigned int wrong = 0;
	//unsigned int main_height = main.get_height();
	//unsigned int main_width = main.get_width();
	unsigned int main_height = main.height;
	unsigned int main_width = main.width;
	unsigned int height = t.height;
	unsigned int width = t.width;
	//int maxwrong = (int) (maxwrong_ratio * t.get_height() * t.get_width()); 
	int maxwrong2 = (int) (maxwrong_ratio * height * width);
	//The next cases are not possible
	PixelStr pixel = main.pixel_data[main_width*h + w];
	PixelStr first = t.pixel_data[0];
	//The next cases are not possible
	/*if(abs(pixel.r - first.r) > maxdif || 
		abs(pixel.g - first.g) > maxdif || 
		abs(pixel.b - first.b) > maxdif)
		return false;*/
	Point first_edge = t.reference;
	
#if 0
	if(abs(pixel.r - first_edge.r) > maxdif || 
		abs(pixel.g - first_edge.g) > maxdif || 
		abs(pixel.b - first_edge.b) > maxdif)
		return false;

	//cout<<w<<" "<<h<<endl;
	//return true;
#endif		
	// check if the current point is also an edge in the main image	
	bool found = false;
	//if(t.min_type=='r')		
		for(Point c : main.edgeR[t.index_min])
			if(match_point(first_edge, c))
				found = true;
	//else if(t.min_type=='g')		
		for(Point c : main.edgeG[t.index_min])
			if(match_point(first_edge, c))
				found = true;
	//else if(t.min_type=='b')		
		for(Point c : main.edgeB[t.index_min])
			if(match_point(first_edge, c))
				found = true;
	if(found == false) return false;
	
	
	float PI = 3.14159265;
	
	{		
		int maxwrong = (int) (maxwrong_ratio_edge * t.rarest->size()); 
		if(maxwrong==0) maxwrong = 1;
		register unsigned int misses = 0;
		int maxwrong3 = (int) (maxwrong_ratio * t.biggest->size()); 
		if(maxwrong3==0) maxwrong3 = 1;
		
		unsigned int iterations=0;
		
		unsigned int size = t.rarest->size();
		unsigned int size2 = t.biggest->size();
		
		PixelStr pixel2;
		PixelStr *pixel3;
		PixelStr *p;
		//Point e;
		unsigned int scale;
		unsigned int maxscale = 1;
		
		unsigned int x, y;
		float angle;		

		int h2, w2;
		float alfa, beta, sx, sy, theta;
		for(alfa = 0.2; alfa <= max_scale; alfa += 0.1)
			//for(beta = 0.2; beta <= max_scale; beta += 0.2)
			{ beta = alfa;
				for(angle = -15; angle <= 345; angle += 1, iterations++) 
					/*for(angle = 0; angle <= 360; angle += 5) 
					for(sx = -0.01; sx <= 0.01; sx += 0.01)
						for(sy = -0.01; sy <= 0.01; sy += 0.01, iterations++)*/						
					{
						//theta = 0;
						sx = 0;
						sy = 0;
						theta = angle * PI / 180.0;
						if(angle >= 15) angle += 4;

						int a11, a12, a21, a22;	
						int factor = 10000;						
						
						misses = 0;		
							
						// apply the transform to the min-edge template points and calculate number of misses 

						//Point *e = &t.edgeR[t.index_min][0];
#if 1
						Point e;
						for(unsigned int k=0; k<size && misses<=maxwrong; k++)	//{
						//for(Point e : t.edgeR[t.index_min])	
						{
							//e = t.edgeR[t.index_min][k];
							e = (*t.rarest)[k];
							Position pos = t.edgetransform[k + iterations*size];
							//x = e.x;
							//y = e.y;								
							int x2, y2; 
							
							x2 = floor(pos.x + h*factor) / factor;
							y2 = floor(pos.y + w*factor) / factor;
							
							if(x2<0 || y2<0 || x2>=main_height || y2>=main_width) misses++;
							else
							{
								//PixelStr pixel2;
								pixel2 = main.pixel_data[main_width*x2 + y2];//pixel2;
								
								if(abs(e.r - pixel2.r) > maxdif || abs(e.g - pixel2.g) > maxdif || abs(e.b - pixel2.b) > maxdif)
									misses++;
							}
						}	
									
						//continue with the next iteration if number of misses is not acceptable							
						if(misses > maxwrong) continue; 
						
						// we now have a candidate match

						a11 =  static_cast<int>(floor((alfa*cos(theta) - alfa*sx*sin(theta)) * factor + 0.5));
						a12 =  static_cast<int>(floor((alfa*sin(theta) + alfa*sx*cos(theta)) * factor + 0.5));
						a21 =  static_cast<int>(floor((beta*sy*cos(theta) - beta*sin(theta)) * factor + 0.5));
						a22 =  static_cast<int>(floor((beta*sy*sin(theta) + beta*cos(theta)) * factor + 0.5));
						
						// apply the transform to the max-edge template points and calculate number of misses
						misses = 0;
						 
						Point e2;
						for(unsigned int k=0; k<size2 && misses<=maxwrong3; k++)	//{
						//for(Point e : t.edgeR[t.index_min])	
						{
							//e = t.edgeR[t.index_min][k];
							e2 = (*t.biggest)[k];
							//Position pos = t.edgetransform[k + iterations*size];
							x = e2.x;
							y = e2.y;								
							int x2, y2; 
							
							x2 = static_cast<int>(floor(factor*h + (a11*(x-first_edge.x) + a12*(y-first_edge.y)))) / factor ;
							y2 = static_cast<int>(floor(factor*w + (a21*(x-first_edge.x) + a22*(y-first_edge.y)))) / factor ;

							if(x2<0 || y2<0 || x2>=main_height || y2>=main_width) misses++;
							else
							{
								//PixelStr pixel2;
								pixel2 = main.pixel_data[main_width*x2 + y2];//pixel2;
								
								if(abs(e2.r - pixel2.r) > maxdif || abs(e2.g - pixel2.g) > maxdif || abs(e2.b - pixel2.b) > maxdif)
									misses++;
							}
						}	
									
						//continue with the next iteration if number of misses is not acceptable							
						if(misses > maxwrong3) continue; 
						
						// end of second filter

						// last filter, compare all template points against their transformed positions in the main image
						misses = 0;
#endif
						misses = 0;

						p = &t.pixel_data[0];						
					
						for(unsigned int x=0; x<height && misses<=maxwrong2; x++)
						{
							for(unsigned int y=0; y<width && misses<=maxwrong2; y++,p++)
							{	
								int x2, y2;
																
								x2 = static_cast<int>(floor(factor*h + (a11*(x-first_edge.x) + a12*(y-first_edge.y)))) / factor ;
								y2 = static_cast<int>(floor(factor*w + (a21*(x-first_edge.x) + a22*(y-first_edge.y)))) / factor ;
								
								if(x2<0 || y2<0 || x2>=main_height || y2>=main_width) misses++;
								else
								{
									pixel3 = &main.pixel_data[main_width*x2 + y2];
								
									if(abs(p->r - pixel3->r) > maxdif || abs(p->g - pixel3->g) > maxdif || abs(p->b - pixel3->b) > maxdif)
										misses++;
								}							
							}
						}
						//cout<<hits<<" hits, max_allowed_missed = "<<maxwrong2<<endl;
						if(misses > maxwrong2) continue;
						
						// we now have a match, but we still need to compute the (0,0) - transform

						h2 = static_cast<int>(floor(factor*h - a11*(first_edge.x) - a12*(first_edge.y) + 0.5)) / factor;
						w2 = static_cast<int>(floor(factor*w - a21*(first_edge.x) - a22*(first_edge.y) + 0.5)) / factor;

						// h and w are passed by reference, thus these values will be also available in the main function
						h = h2;
						w = w2;
						//cout<<alfa<<" "<<beta<<" "<<angle<<" "<<sx<<" "<<sy<<" "<<h<<" "<<w<<endl;
						return true;
							//iterations++;
					}
					
		//cout<<iterations<<endl;
			}// end of alfa-beta loop
	}
	// if we get here, it means no matches were found 
	return false;
}



bool match_template(const Accelerate::Image& main_image, const Accelerate::Image& template_image, unsigned int &h, unsigned int &w, int max_scale){
#if 0
	unsigned int minh, maxh, minw, maxw;
	minh = 0;//400;
	maxh = 1000;//800;

	minw = 0;//3200;	
	maxw = 1000;//3800;
	if(h>maxh || w>maxw || h<minh || w<minw) return false;
#endif
	return match1(main_image, template_image, h, w, max_scale);
	//return false;
}

bool match_point(const Point& p1, const Point& p2){	
	if(abs(p1.r - p2.r) > maxdif) return false;
	if(abs(p1.g - p2.g) > maxdif) return false;
	if(abs(p1.b - p2.b) > maxdif) return false;
	return true;
}




