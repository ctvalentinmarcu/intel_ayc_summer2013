/*!
 * \file BitmapImporter.cpp
 * \brief Contains the implementation of the functions declared in BitmapImporter.h
 * \author Cédric Andreolli (Intel)
 * \date 10 April 2013
 */
#include "BitmapImporter.h" 

using namespace std;
using namespace Accelerate;

Accelerate::Image Accelerate::Image::create_image_from_bitmap(const std::string file_name, int max_scale, unsigned int main_height, unsigned int main_width){

	unsigned char maxdif = 10;
	float maxwrong_ratio = 0.01;

	Accelerate::Image result;
	ifstream file(file_name.c_str(), ios::in|ios::binary|ios::ate);
	if (file.is_open())
  	{
		Accelerate::HeaderStr header;
    		file.seekg (0, ios::beg);
		//Read the header first
    		file.read ((char*)&header, sizeof(HeaderStr));

		result.width = header.width;
		result.height = header.height;

		//Put the cursor on the BMP data
		file.seekg(header.offset, ios::beg);
		int bytes_per_pixel = header.bit_per_pixel / 8;
		int padding = ((header.width*bytes_per_pixel) % 4 == 0) ? 0 : 4 - ((header.width*bytes_per_pixel) % 4);
		//Allocate the size needed to read the BMP data
		int size = sizeof(char)*(header.height*(header.width + padding))*bytes_per_pixel;
		unsigned char* data = (unsigned char*) malloc(size);
		//Read the data
		file.read((char*)data, size);
		//Create the Accelerate::Bitmap object
		result.pixel_data = (PixelStr*) malloc(header.width*header.height*sizeof(PixelStr));
		result.aux = (PixelStrInt*) malloc(header.width*header.height*sizeof(PixelStrInt));
		int hr[256], hg[256], hb[256]; // histogram
		unsigned int offset = 0;
		//In the Bitmap format, pixels are in a reversed order
		for(int i=header.height-1; i>=0; i--){
			for(int j=0; j<header.width; j++){
				unsigned char b = data[offset++];
				unsigned char g = data[offset++];
				unsigned char r = data[offset++];
				/*int b = data[offset++];
				int g = data[offset++];
				int r = data[offset++];*/
				result.pixel_data[i*header.width + j].b = b;//(unsigned char) (b / maxdif * maxdif);
				result.pixel_data[i*header.width + j].g = g;//(unsigned char) (g / maxdif * maxdif);
				result.pixel_data[i*header.width + j].r = r;//(unsigned char) (r / maxdif * maxdif);
			}
			offset+=padding;			
		}

		result.maxwrong = (int) (maxwrong_ratio * header.height * header.width);

		bool isEdgeR=false, isEdgeG=false, isEdgeB=false;

#if 1
		//result.edgeR = new vector<vector<int>>();
		for(int i=0; i<256; i++){
			//vector<int> v();
			result.edgeR.push_back(vector<Point>());
			result.edgeG.push_back(vector<Point>());
			result.edgeB.push_back(vector<Point>());
		}

#if 0
		// apply normalization
		int minr, maxr, ming, maxg, minb, maxb;
		minr = ming = minb = 256;
		maxr = maxg = maxb = -1;
		for(int i=0; i<256 && minr<0; i++)
			if(hr[i]>0) minr = i;
		for(int i=0; i<256 && ming<0; i++)
			if(hg[i]>0) ming = i;
		for(int i=0; i<256 && minb<0; i++)
			if(hb[i]>0) minb = i;

		for(int i=255; i>=0 && maxr>255; i--)
			if(hr[i]>0) maxr = i;
		for(int i=255; i>=0 && maxg>255; i--)
			if(hg[i]>0) maxg = i;
		for(int i=255; i>=0 && maxb>255; i--)
			if(hb[i]>0) maxb = i;
		
		
		for(int i=0; i<header.height*header.width; i++){
			float t;
			t = ((float) (result.pixel_data[i].r - minr)) * 255 / (maxr - minr);
			result.pixel_data[i].r = (unsigned char) floor(t + 0.5);
			t = ((float) (result.pixel_data[i].g - ming)) * 255 / (maxg - ming);
			result.pixel_data[i].g = (unsigned char) floor(t + 0.5);
			t = ((float) (result.pixel_data[i].b - minb)) * 255 / (maxb - minb);
			result.pixel_data[i].b = (unsigned char) floor(t + 0.5);
		}
#endif


		// apply histogram equalization
#if 0		
		for(int i=1; i<256; i++){
			hr[i] += hr[i-1];
			hg[i] += hg[i-1];
			hb[i] += hb[i-1];
		}
		
		for(int i=0; i<header.height*header.width; i++){
			float t;
			t = ((float) hr[result.pixel_data[i].r] - hr[0]) * 255 / (header.height*header.width - hr[0]); 
			result.pixel_data[i].r = (unsigned char) floor(t + 0.5);
			t = ((float) hg[result.pixel_data[i].g] - hg[0]) * 255 / (header.height*header.width - hg[0]); 
			result.pixel_data[i].g = (unsigned char) floor(t + 0.5);
			t = ((float) hb[result.pixel_data[i].b] - hb[0]) * 255 / (header.height*header.width - hb[0]); 
			result.pixel_data[i].b = (unsigned char) floor(t + 0.5);
		}
#endif	


#if 0
		// apply the gaussian filter using a 5x5 kernel
		for(int i=0; i<header.height; i++){  
			for(int j=0; j<header.width; j++){
				result.aux[i*header.width + j].r = result.pixel_data[i*header.width + j].r;
				result.aux[i*header.width + j].g = result.pixel_data[i*header.width + j].g;
				result.aux[i*header.width + j].b = result.pixel_data[i*header.width + j].b;
				if(i>=2 && i<=header.height-3 && j>=2 && j<=header.width-3){
					float auxr, auxg, auxb;
					auxr = result.pixel_data[(i-2)*header.width + j-2].r + 4*result.pixel_data[(i-2)*header.width + j-1].r + 
							7*result.pixel_data[(i-2)*header.width + j].r + 4*result.pixel_data[(i-2)*header.width + j+1].r + result.pixel_data[(i-2)*header.width + j+2].r;
					auxr += 4*result.pixel_data[(i-1)*header.width + j-2].r + 16*result.pixel_data[(i-1)*header.width + j-1].r + 
							26*result.pixel_data[(i-1)*header.width + j].r + 16*result.pixel_data[(i-1)*header.width + j+1].r + 4*result.pixel_data[(i-1)*header.width + j+2].r;
					auxr += 7*result.pixel_data[(i)*header.width + j-2].r + 26*result.pixel_data[(i)*header.width + j-1].r + 
							41*result.pixel_data[(i)*header.width + j].r + 26*result.pixel_data[(i)*header.width + j+1].r + 7*result.pixel_data[(i)*header.width + j+2].r;
					auxr += 4*result.pixel_data[(i+1)*header.width + j-2].r + 16*result.pixel_data[(i+1)*header.width + j-1].r + 
							26*result.pixel_data[(i+1)*header.width + j].r + 16*result.pixel_data[(i+1)*header.width + j+1].r + 4*result.pixel_data[(i+1)*header.width + j+2].r;
					auxr += result.pixel_data[(i+2)*header.width + j-2].r + 4*result.pixel_data[(i+2)*header.width + j-1].r + 
							7*result.pixel_data[(i+2)*header.width + j].r + 4*result.pixel_data[(i+2)*header.width + j+1].r + result.pixel_data[(i+2)*header.width + j+2].r;

					auxg = result.pixel_data[(i-2)*header.width + j-2].g + 4*result.pixel_data[(i-2)*header.width + j-1].g + 
							7*result.pixel_data[(i-2)*header.width + j].g + 4*result.pixel_data[(i-2)*header.width + j+1].g + result.pixel_data[(i-2)*header.width + j+2].g;
					auxg += 4*result.pixel_data[(i-1)*header.width + j-2].g + 16*result.pixel_data[(i-1)*header.width + j-1].g + 
							26*result.pixel_data[(i-1)*header.width + j].g + 16*result.pixel_data[(i-1)*header.width + j+1].g + 4*result.pixel_data[(i-1)*header.width + j+2].g;
					auxg += 7*result.pixel_data[(i)*header.width + j-2].g + 26*result.pixel_data[(i)*header.width + j-1].g + 
							41*result.pixel_data[(i)*header.width + j].g + 26*result.pixel_data[(i)*header.width + j+1].g + 7*result.pixel_data[(i)*header.width + j+2].g;
					auxg += 4*result.pixel_data[(i+1)*header.width + j-2].g + 16*result.pixel_data[(i+1)*header.width + j-1].g + 
							26*result.pixel_data[(i+1)*header.width + j].g + 16*result.pixel_data[(i+1)*header.width + j+1].g + 4*result.pixel_data[(i+1)*header.width + j+2].g;
					auxg += result.pixel_data[(i+2)*header.width + j-2].g + 4*result.pixel_data[(i+2)*header.width + j-1].g + 
							7*result.pixel_data[(i+2)*header.width + j].g + 4*result.pixel_data[(i+2)*header.width + j+1].g + result.pixel_data[(i+2)*header.width + j+2].g;

					auxb = result.pixel_data[(i-2)*header.width + j-2].b + 4*result.pixel_data[(i-2)*header.width + j-1].b + 
							7*result.pixel_data[(i-2)*header.width + j].b + 4*result.pixel_data[(i-2)*header.width + j+1].b + result.pixel_data[(i-2)*header.width + j+2].b;
					auxb += 4*result.pixel_data[(i-1)*header.width + j-2].b + 16*result.pixel_data[(i-1)*header.width + j-1].b + 
							26*result.pixel_data[(i-1)*header.width + j].b + 16*result.pixel_data[(i-1)*header.width + j+1].b + 4*result.pixel_data[(i-1)*header.width + j+2].b;
					auxb += 7*result.pixel_data[(i)*header.width + j-2].b + 26*result.pixel_data[(i)*header.width + j-1].b + 
							41*result.pixel_data[(i)*header.width + j].b + 26*result.pixel_data[(i)*header.width + j+1].b + 7*result.pixel_data[(i)*header.width + j+2].b;
					auxb += 4*result.pixel_data[(i+1)*header.width + j-2].b + 16*result.pixel_data[(i+1)*header.width + j-1].b + 
							26*result.pixel_data[(i+1)*header.width + j].b + 16*result.pixel_data[(i+1)*header.width + j+1].b + 4*result.pixel_data[(i+1)*header.width + j+2].b;
					auxb += result.pixel_data[(i+2)*header.width + j-2].b + 4*result.pixel_data[(i+2)*header.width + j-1].b + 
							7*result.pixel_data[(i+2)*header.width + j].b + 4*result.pixel_data[(i+2)*header.width + j+1].b + result.pixel_data[(i+2)*header.width + j+2].b;

					result.aux[i*header.width + j].r = (unsigned char) ( auxr / 273);// / maxdif * maxdif;
					result.aux[i*header.width + j].g = (unsigned char) ( auxg / 273);// / maxdif * maxdif;
					result.aux[i*header.width + j].b = (unsigned char) ( auxb / 273);// / maxdif * maxdif;
				}
			}
		}

		for(int i=0; i<header.height; i++)  
			for(int j=0; j<header.width; j++){
				result.pixel_data[i*header.width + j].r = (unsigned char) result.aux[i*header.width + j].r;
				result.pixel_data[i*header.width + j].g = (unsigned char) result.aux[i*header.width + j].g;
				result.pixel_data[i*header.width + j].b = (unsigned char) result.aux[i*header.width + j].b;
			}
		// end of gaussian filter
#endif


#if 0
		// apply the LoG filter using a 3x3 kernel
		for(int i=0; i<header.height; i++){  
			for(int j=0; j<header.width; j++){
				//result.aux[i*header.width + j] = result.pixel_data[i*header.width + j];
				if(i>=1 && i<=header.height-2 && j>=1 && j<=header.width-2){
					int auxr, auxg, auxb;
					
					auxr = -result.pixel_data[(i-1)*header.width + j-1].r - result.pixel_data[(i-1)*header.width + j].r - result.pixel_data[(i-1)*header.width + j+1].r;
					auxr += -result.pixel_data[(i)*header.width + j-1].r + 8*result.pixel_data[(i)*header.width + j].r - result.pixel_data[(i)*header.width + j+1].r;
					auxr += -result.pixel_data[(i+1)*header.width + j-1].r - result.pixel_data[(i+1)*header.width + j].r - result.pixel_data[(i+1)*header.width + j+1].r;					

					auxg = -result.pixel_data[(i-1)*header.width + j-1].g - result.pixel_data[(i-1)*header.width + j].g - result.pixel_data[(i-1)*header.width + j+1].g;
					auxg += -result.pixel_data[(i)*header.width + j-1].g + 8*result.pixel_data[(i)*header.width + j].g - result.pixel_data[(i)*header.width + j+1].g;
					auxg += -result.pixel_data[(i+1)*header.width + j-1].g - result.pixel_data[(i+1)*header.width + j].g - result.pixel_data[(i+1)*header.width + j+1].g;

					auxb = -result.pixel_data[(i-1)*header.width + j-1].b - result.pixel_data[(i-1)*header.width + j].b - result.pixel_data[(i-1)*header.width + j+1].b;
					auxb += -result.pixel_data[(i)*header.width + j-1].b + 8*result.pixel_data[(i)*header.width + j].b - result.pixel_data[(i)*header.width + j+1].b;
					auxb += -result.pixel_data[(i+1)*header.width + j-1].b - result.pixel_data[(i+1)*header.width + j].b - result.pixel_data[(i+1)*header.width + j+1].b;
					/*
					if(auxr > maxdif) isEdgeR = true;
					if(auxg > maxdif) isEdgeG = true;
					if(auxb > maxdif) isEdgeB = true;
					*/
					/*auxr += result.pixel_data[i*header.width + j].r;
					auxg += result.pixel_data[i*header.width + j].g;
					auxb += result.pixel_data[i*header.width + j].b;*/
					
					/*if(auxr < 0) auxr = 0;
					if(auxr > 255) auxr = 255;
					if(auxg < 0) auxg = 0;
					if(auxg > 255) auxg = 255;
					if(auxb < 0) auxb = 0;
					if(auxb > 255) auxb = 255;*/
					result.aux[i*header.width + j].r =  ( auxr );// / maxdif * maxdif;
					result.aux[i*header.width + j].g =  ( auxg );// / maxdif * maxdif;
					result.aux[i*header.width + j].b =  ( auxb );// / maxdif * maxdif;
					
				}
			}
		}
		
		// end of LoG filter
#endif

		Point *p = new Point();
		PixelStr n1, n2, n3, n4, n5, n6, n7, n8;
		PixelStr *pixel;
		//pixel = (PixelStr*) malloc(sizeof(PixelStr));
		pixel = new PixelStr();
		// group points by colour
		for(int i=0; i<header.height; i++){  //cout<<endl;
			for(int j=0; j<header.width; j++){
				//PixelStr *pixel;
				//pixel = (PixelStr*) malloc(sizeof(PixelStr));
				pixel->r = result.pixel_data[i*header.width + j].r;
				pixel->g = result.pixel_data[i*header.width + j].g;
				pixel->b = result.pixel_data[i*header.width + j].b;					
				unsigned int nr = 0;
				//unsigned int pixelr = pixel->r;
				//cout<<pixelr<<"   ";

				
				bool isEdgeR=false, isEdgeG=false, isEdgeB=false;
				
				// the 8 neighbors
				//PixelStr n1, n2, n3, n4, n5, n6, n7, n8;
				if(i>0){
					if(j>0)
						n1 = result.pixel_data[(i-1)*header.width + j-1];
						n2 = result.pixel_data[(i-1)*header.width + j];
					if(j<header.width-1)
						n3 = result.pixel_data[(i-1)*header.width + j+1];
				}
				if(i<header.height-1) {
					if(j>0)
						n4 = result.pixel_data[(i+1)*header.width + j-1];
						n5 = result.pixel_data[(i+1)*header.width + j];
					if(j<header.width-1)
						n6 = result.pixel_data[(i+1)*header.width + j+1];
				}
				if(j>0) 
					n7 = result.pixel_data[(i)*header.width + j-1];
				if(j<header.width-1) 
					n8 = result.pixel_data[(i)*header.width + j+1];
				
				//isEdgeR = isEdgeG = isEdgeB = false;
				
				isEdgeR = (i>0) && ((j>0 && abs(n1.r - pixel->r)>= maxdif) || (abs(n2.r - pixel->r)>= maxdif) || (j<header.width-1 && abs(n3.r - pixel->r)>= maxdif) ); 
				isEdgeR = isEdgeR || ( (i<header.height-1) && ( (j>0 && abs(n4.r - pixel->r)>=maxdif)||(abs(n5.r - pixel->r)>=maxdif)||(j<header.width-1 && abs(n6.r - pixel->r)>=maxdif) ) );
				isEdgeR = isEdgeR || ( (j>0 && abs(n7.r - pixel->r) >= maxdif) || (j<header.width-1 && abs(n8.r - pixel->r) >= maxdif) );

				isEdgeG = (i>0) && ((j>0 && abs(n1.g - pixel->g)>= maxdif) || (abs(n2.g - pixel->g)>= maxdif) || (j<header.width-1 && abs(n3.g - pixel->g)>= maxdif) ); 
				isEdgeG = isEdgeG || ( (i<header.height-1) && ( (j>0 && abs(n4.g - pixel->g)>=maxdif)||(abs(n5.g - pixel->g)>=maxdif)||(j<header.width-1 && abs(n6.g - pixel->g)>=maxdif) ) );
				isEdgeG = isEdgeG || ( (j>0 && abs(n7.g - pixel->g) >= maxdif) || (j<header.width-1 && abs(n8.g - pixel->g) >= maxdif) );

				isEdgeB = (i>0) && ((j>0 && abs(n1.b - pixel->b)>= maxdif) || (abs(n2.b - pixel->b)>= maxdif) || (j<header.width-1 && abs(n3.b - pixel->b)>= maxdif) ); 
				isEdgeB = isEdgeB || ( (i<header.height-1) && ( (j>0 && abs(n4.b - pixel->b)>=maxdif)||(abs(n5.b - pixel->b)>=maxdif)||(j<header.width-1 && abs(n6.b - pixel->b)>=maxdif) ) );
				isEdgeB = isEdgeB || ( (j>0 && abs(n7.b - pixel->b) >= maxdif) || (j<header.width-1 && abs(n8.b - pixel->b) >= maxdif) );
				
#if 0
				int minr, ming, minb, maxr, maxg, maxb;
				int treshold = 0;
				//for(int i=0; i<header.height; i++){  
					//for(int j=0; j<header.width; j++){
						//result.aux[i*header.width + j] = result.pixel_data[i*header.width + j];
						if(i>=1 && i<=header.height-2 && j>=1 && j<=header.width-2){
							minr = ming = minb = 255;
							maxr = maxg = maxb = -255;
							for(int ii = i-1; ii <= i+1; ii++)
								for(int jj = j-1; jj <= j+1; jj++){
									minr = (minr < result.aux[ii*header.width + jj].r) ? minr : result.aux[ii*header.width + jj].r;
									ming = (ming < result.aux[ii*header.width + jj].g) ? ming : result.aux[ii*header.width + jj].g;
									minb = (minb < result.aux[ii*header.width + jj].b) ? minb : result.aux[ii*header.width + jj].b;

									maxr = (maxr > result.aux[ii*header.width + jj].r) ? maxr : result.aux[ii*header.width + jj].r;
									maxg = (maxg > result.aux[ii*header.width + jj].g) ? maxg : result.aux[ii*header.width + jj].g;
									maxb = (maxb > result.aux[ii*header.width + jj].b) ? maxb : result.aux[ii*header.width + jj].b;
								}
							if(maxr > 0 && minr < 0 && maxr-minr > treshold) isEdgeR = true;
							if(maxg > 0 && ming < 0 && maxg-ming > treshold) isEdgeG = true;
							if(maxb > 0 && minb < 0 && maxb-minb > treshold) isEdgeB = true;
						}
				//	}
				//}
#endif
				int index = i*header.width+j;
				int indexr, indexg, indexb;

				// approximate the intensity by a factor of maxdif
				// indexr...indexg will only have intensity values of 0, maxdif, 2*maxdif etc
				// but the pixel_data array will still contain the original values
				indexr = (int) (result.pixel_data[index].r) / maxdif * maxdif;
				indexg = (int) (result.pixel_data[index].g) / maxdif * maxdif;
				indexb = (int) (result.pixel_data[index].b) / maxdif * maxdif;

				/*indexr = (int) (result.pixel_data[index].r);
				indexg = (int) (result.pixel_data[index].g);
				indexb = (int) (result.pixel_data[index].b);*/
			
				// construct the edge point
				p->x = i;
				p->y = j;
				p->r = pixel->r;
				p->g = pixel->g;
				p->b = pixel->b;
				
				// put the edge point into the corresponding vector(s)
				if(isEdgeR==true) result.edgeR[indexr].push_back(*p);
				if(isEdgeG==true) result.edgeG[indexg].push_back(*p);
				if(isEdgeB==true) result.edgeB[indexb].push_back(*p);
				
			}			
		}

		//find the rarest coloured edge
		int min = 10000;

		for(int i=0; i<256; i++){
			int sizer = result.edgeR[i].size();
			int sizeg = result.edgeG[i].size();
			int sizeb = result.edgeB[i].size();
			if(result.edgeR[i].size()>1 && result.edgeR[i].size() < min) //if(abs(i - index_max_hist_r) > maxdif)
			{
				min = result.edgeR[i].size();
				result.index_min = i;
				result.min_type = 'r';
				result.rarest = &result.edgeR[i];
			}
			if(result.edgeG[i].size()>1 && result.edgeG[i].size() < min) //if(abs(i - index_max_hist_g) > maxdif)
			{
				min = result.edgeG[i].size();
				result.index_min = i;
				result.min_type = 'g';
				result.rarest = &result.edgeG[i];
			}
			if(result.edgeB[i].size()>1 && result.edgeB[i].size() < min) //if(abs(i - index_max_hist_b) > maxdif)
			{
				min = result.edgeB[i].size();
				result.index_min = i;
				result.min_type = 'b';
				result.rarest = &result.edgeB[i];
			}
		}

		// find the "largest" monocolor edge
		int max = 0;
		
		for(int i=0; i<256; i++){
			if(result.edgeR[i].size()>1 && result.edgeR[i].size() > max){
				max = result.edgeR[i].size();
				result.index_max = i;
				result.max_type = 'r';
				result.biggest = &result.edgeR[i];
			}
			if(result.edgeG[i].size()>1 && result.edgeG[i].size() > max){
				max = result.edgeG[i].size();
				result.index_max = i;
				result.max_type = 'g';
				result.biggest = &result.edgeG[i];
			}
			if(result.edgeB[i].size()>1 && result.edgeB[i].size() > max){
				max = result.edgeB[i].size();
				result.index_max = i;
				result.max_type = 'b';
				result.biggest = &result.edgeB[i];
			}
		}

		result.reference = (*result.rarest)[0];

		//cout<<result.index_min<<result.min_type<<min<<"; reference: "<<result.reference.x<<" "<<result.reference.y<<endl;
		
		float angle;
		float theta;
		float alfa, beta, sx, sy;
		int a11, a12, a21, a22;
		float PI = 3.14159265;
		
		int a, b;
		//int max_scale = 4;
		int size2 = result.rarest->size();
		Point e;
		int x, y;
		//result.edgetransform = vector<Position>();
		Position *pos = new Position();

		//max_scale = 1;
		for(alfa = 0.2; alfa <= max_scale; alfa += 0.1)
			//for(beta = 0.2; beta <= max_scale; beta += 0.2)
				for(angle = -15; angle <= 345; angle += 1) 
					//for(sx = -0.01; sx <= 0.01; sx += 0.01)
						//for(sy = -0.01; sy <= 0.01; sy += 0.01)
				{
					sx = 0;
					sy = 0;
					beta = alfa;
					theta = angle * PI / 180.0;
					if(angle >= 15) angle += 4;

					int factor = 10000;
					a11 =  static_cast<int>(floor((alfa*cos(theta) - alfa*sx*sin(theta)) * factor + 0.5));
					a12 =  static_cast<int>(floor((alfa*sin(theta) + alfa*sx*cos(theta)) * factor + 0.5));
					a21 =  static_cast<int>(floor((beta*sy*cos(theta) - beta*sin(theta)) * factor + 0.5));
					a22 =  static_cast<int>(floor((beta*sy*sin(theta) + beta*cos(theta)) * factor + 0.5));
					
					for(unsigned int k=0; k<size2; k++){
						e = (*result.rarest)[k];
						x = e.x;
						y = e.y;	
						//cout<<x<<" "<<y<<endl;
						int x2, y2; 
						
						x2 = a11*(x-result.reference.x) + a12*(y-result.reference.y);
						y2 = a21*(x-result.reference.x) + a22*(y-result.reference.y);

						pos->x = x2;
						pos->y = y2;
						result.edgetransform.push_back(*pos);
						
					}
					
				}
		
		
#endif
		file.close();
		free(data);
	}
	return result;
}

Accelerate::Image::Image(){
		
}

Accelerate::Image::~Image(){
	free(pixel_data);
	//free(edgeR);
	//free(edgeG);
	//free(edgeB);
}

Accelerate::PixelStr Accelerate::Image::get_pixel(unsigned int i, unsigned int j) const{
	return this->pixel_data[this->width*i + j];
}

std::ostream& operator<<(std::ostream &o, const Accelerate::PixelStr& p){
    return o << "[" << (int)p.r << ", " << (int)p.g << ", " << (int)p.b  << "] ";
}

Accelerate::Image Accelerate::Image::scale_image(unsigned int scale) const{
	Accelerate::Image result;

	result.width = width * scale;
	result.height = height * scale;
	
	result.pixel_data = (PixelStr*) malloc(result.width*result.height*sizeof(PixelStr));

	for(unsigned int w=0; w<this->width; w++){
		for(unsigned int i=0; i<scale; i++)
			copy_column(result, w, scale); 
	}
	
	return result;
}

void Accelerate::Image::copy_column(Accelerate::Image& result, unsigned int column_number, unsigned int scale) const{
	//retrieve the column indice in the result image
	unsigned int first_column_indice = column_number * scale;
	for(unsigned int i=0; i<scale; i++){
		for(unsigned int row=0; row<height; row++){
			for(unsigned int j=0; j<scale; j++){
				result.pixel_data[(row*scale + j)*result.width + first_column_indice + i] = this->pixel_data[row*this->width + column_number]; 
			}
		}
	}
} 


std::ostream& operator<<(std::ostream &o, Accelerate::Image& im){
	Accelerate::PixelStr* pixels = im.get_pixels();	
	for(unsigned int i=0; i<im.get_height(); i++){
		for(unsigned int j=0; j<im.get_width(); j++){
			Accelerate::PixelStr pixel = pixels[i*im.get_width() + j];
			o<<pixel<<" ";
		}
		o<<std::endl;
	}


	return o; 
}
