#ifndef PMT_h
#define PMT_h

#include <vector>

#include "Pixel.h"
#include "PixelSpot.h"

class PMT
{

public:
                PMT(int id){Init(id);};
                PMT(){Init(-1);};
                ~PMT();

                int GetId(){return id;};

                Pixel * GetPixel(int id,int create=1);

                std::vector<int>  GetPixelVector();

                void Dump();
                
                int GetNStrips(){return nstrips;};
                void IncStripCount(){nstrips++;};
                
                int GetStatus(){return status;};
                void SetStatus(int s){status=s;};
                
                static int printit;

private:
                void Init(int id);
        
                int id;
                int nstrips;
                
                std::vector<Pixel> pixels;
                int status;
};

#endif
