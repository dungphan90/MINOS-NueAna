#ifndef PixelSpot_h
#define PixelSpot_h

#include <vector>

class Pixel;

class PixelSpot
{
        public:
                PixelSpot();
                PixelSpot(int id);
                ~PixelSpot();

                int GetId(){return id;};

                void AddStrip(int strip, double q, double strip_e);
                
                void Dump();
                
                void SetPixel(Pixel * p);
                Pixel * GetPixel();
                
                double GetE(){return sum_q;};
                double GetStripE(){return sum_strip_e;};
                int GetNStrips(){return n;};
                
                std::vector<int> GetStrips(){return strips;};
                double GetStripQ(int idx){if(idx>-1&&idx<n)return q[idx];return 0;};
                
                void AddE(double e);
                
                void TakeE(double e);
                
                int GetStatus(){return status;};

                static int printit;
                
        private:
                void Init(int id);      
                int id;
                std::vector<int> strips;
                std::vector<double> q;
                std::vector<double>strip_e;
                double sum_q;
                double sum_strip_e;
                double spot_e;
                int n;
                Pixel * parentPixel;
// keep the status of the pixelspot... if its 0, we didn't do anything to it!
                int status;
};

#endif
