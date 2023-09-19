#include "MFileReader.h"
#include "StringTokenizer.h"
#include "SolidDelegate.h"
#include <list>		//list of std

using namespace MeshLib;


MFileReader::MFileReader ()
{
	vID = 0;
	fID = 0;
};
		

MFileReader::~MFileReader()
{

};
		
//!	read input file, handle different contents
void MFileReader::readToSolid( Solid *mesh, std::istream &in)
{
	SolidDelegate delegate;
	char strLine [1024]  = {0};
	
	while (in && !in.eof() && in.getline(strLine, 1024)) 
	{

		if( strlen( strLine ) == 0 ) continue;

		StringTokenizer tokenizer( strLine, " \r\n\t" );

		List<char> & tokens = tokenizer.tokens();
		ListNode<char> * node = tokens.head();

		char * token = node->data();

		//!read vertex coordinates, generate vertex
		if( strcmp( token, "Vertex" ) == 0 )
		{
			Point p;
			int v_id;
			node= node->next();
			if( strlen(node->data()) != 0)
			{
				v_id = atoi( node->data() );
			}

			//read coordinates
			for(int i=0; i<3; i++)
			{
				node = node->next();
				if( strlen(node->data()) != 0)
				{
					p[i] = atof( node->data() );
				}
				else i--;
			}
		
			//generate vertex
			Vertex * v = delegate.createVertex(mesh,  v_id);
			v->point() = p;
			v->id()    = v_id;
		}
		
		//!read face information, generate faces
		else if ( strcmp( token, "Face" ) == 0 )
		{
			char * subToken;
			int id =0, f_id=0;
			std::list<int> polygonVertex;

			node =node ->next();
			if( strlen(node->data()) != 0)
			{
				f_id = atoi( node->data() );
			}
			node =node ->next();

			while( node != NULL )
			{
				if( strlen(node->data()) != 0)	//fill empty node
				{
					id = atoi( node ->data() );
					polygonVertex.push_back(id);
				}
				node = node->next();
			}

			int ids[3] = {0};
			ids[0] = polygonVertex.front();
			polygonVertex.pop_front();
			ids[2] = polygonVertex.front();
			polygonVertex.pop_front();

			while(!polygonVertex.empty())
			{
				ids[1] = ids[2];
				ids[2] = polygonVertex.front();
				polygonVertex.pop_front();
				delegate.createFace(mesh, ids, f_id);

			}
		}

		else
		{
		}
	}//end of while

	mesh->labelBoundaryEdges();
	mesh->removeDanglingVertices();
};

//how about w value, the format?
////////////////////////////////
//...............................

//! double to string convertor
std::string 
MFileReader::d2String( double value)
{
	char buffer [20];
	sprintf(buffer, "%g", value);
//	_gcvt (value,6,buffer);
	return buffer;
};

