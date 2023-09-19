#ifndef _MFILEREADER_H
#define _MFILEREADER_H


#include "Solid.h"
#include "Array.h"
#include "Point.h"


namespace MeshLib
{
	//!	 MFileReader.
	/*!	This class defines reading file in format m/n
		and generating corresponding mesh file/n
	*/

	class MFileReader
	{
	private:
		//!	vertex id
		int vID;
		//!	face id
		int fID;

	public:
		//!	 Constructor
		MFileReader ();
		
		//!	Destructor
		~MFileReader();

		//!	read input file, handle different contents
		void readToSolid( Solid *mesh, std::istream &in);

		//how about w value, the format?
		////////////////////////////////
		//...............................

		//! double to string convertor
		std::string d2String( double value);

	};
}
#endif

