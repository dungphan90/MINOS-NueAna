#ifndef Pixel_h
#define Pixel_h

#include <vector>

class PixelSpot;
class PMT;

class Pixel
{
        public:
                Pixel();
                Pixel(int id);
                ~Pixel();

                int GetId(){return id;};

                PixelSpot * GetPixelSpot(int id, int create=1);

                std::vector<int>  GetPixelSpotVector();
                        
                double GetE(){return e;};
                double GetStripE(){return stripE;};
                int GetNStrips();
                
                void AddE(double ee){e+=ee; if(e<0)e=0;};
                void TakeE(double ee);
                        
                void AddStripE(double ee){stripE+=ee;if(stripE<0)stripE=0;};
                
        
                void Dump();

                void SetPMT(PMT* pmt){parentPMT=pmt;};
                PMT* GetPMT(){return parentPMT;};

                int GetStatus(){return status;};
                void SetStatus(int s){status=s;};
                
                int type;
                int hopstomixed;
                
                static int printit;

        private:
                int id;
                std::vector<PixelSpot> pixelspots;
                double e;
                double stripE;
                PMT* parentPMT;
                
                int status;

};

#endif
