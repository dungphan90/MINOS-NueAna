
 /**
 *
 * $Id: StringUtil.h,v 1.1 2005/06/04 15:47:22 msanchez Exp $
 *
 * \class StringUtil
 *
 * \package NueAna
 *
 * \brief A class to do string manipulation
 *
 * Contact: M. Sanchez
 *
 * Created on: Sat Jun  4 01:44:11 2005
 *
 */

#ifndef STRINGUTIL_H
#define STRINGUTIL_H
 
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

namespace StringUtil
{
inline std::string trim(const std::string& s){
  if(s.length() == 0)
    return s;
  std::size_t beg = s.find_first_not_of(" \a\b\f\n\r\t\v");
  std::size_t end = s.find_last_not_of(" \a\b\f\n\r\t\v");
  if(beg == std::string::npos)
    return "";
  return std::string(s, beg, end - beg + 1);
}

inline static int SplitString(const std::string& input,
                              const std::string& delimiter, 
                              vector<std::string>& results)
{
  int iPos = 0;
  int newPos = -1;
  int sizeS2 = delimiter.size();
  int isize = input.size();

  vector<Int_t> positions;

  newPos = input.find (delimiter, 0);

  if( newPos < 0 ) { return 0; }

  int numFound = 0;

  while( newPos > iPos )
  {
    numFound++;
    positions.push_back(newPos);
    iPos = newPos;
    newPos = input.find (delimiter, iPos+sizeS2+1);
  }

  for( int i=0; i <= (int) positions.size(); i++ )
  {
    string s;
    if( i == 0 ) { s = input.substr( i, positions[i] ); }
    int offset = positions[i-1] + sizeS2;
    if( offset < isize )
    {
      if( i == (int) positions.size() )
      {
        s = input.substr(offset);
      }
      else if( i > 0 )
      {
        s = input.substr( positions[i-1] + sizeS2,
          positions[i] - positions[i-1] - sizeS2 );
      }
    }
    if( s.size() > 0 )
    {
      results.push_back(s);
    }
  }
  return numFound;
}
}
#endif // #ifdef STRINGUTIL_H

