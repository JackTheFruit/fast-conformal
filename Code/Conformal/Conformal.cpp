/*
Author: John Varghese George
Description: This program computes a spherical conformal mapping for
the given mesh input. The program utilizes the fast optimization
technique using spherical constraints given in the paper
"Folding-Free Global Conformal Mapping for Genus-0 Surfaces by Harmonic 
Energy Minimization" by Lai et al. and also incorporates Barzilai–Borwein
step size computation for further speed.
*/
#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "MFileReader.h"
#include "SolidDelegate.h"
#include "Solid.h"
#include "iterators.h"
#include <list>

using namespace MeshLib;

//Computes string constant Kuv of edge
double calKuv(Edge* e){
	Vertex *a,*b,*c,*d;
	e->get_vertices(a,b);
	c = e->halfedge(0)->he_next()->target();	//Get vertices opposite to edge
	d = e->halfedge(1)->he_next()->target();
	Point ca,cb,da,db;
	ca = (a->point())-(c->point());				//Get vectors from c
	cb = (b->point())-(c->point());
	da = (a->point())-(d->point());				//Get vectors from d
	db = (b->point())-(d->point());
	double kuv;
	kuv = (ca*cb)/((ca^cb).norm());				//Compute cot(alpha)
	kuv += (da*db)/((da^db).norm());			//Compute cot(alpha)+cot(beta)
	kuv *= 0.5;									//Compute kuv=0.5*(cot(alpha)+cot(beta))
	return(kuv);
}

//Compute energy of harmonic mapping. If kuvMap is not passed,
//Tuette energy is computed
double calEnergy(SolidEdgeIterator* eiter,std::map<Edge*,double>* kuvMap=NULL){
	double energy = 0.0;
	for((*eiter).reset();!(*eiter).end();++(*eiter)){			//Aggregate energy contribution from each edge
		Edge *e = *(*eiter);
		Vertex *u,*v;
		e->get_vertices(u,v);
		Point fu = u->point();
		Point fv = v->point();
		if(kuvMap)
			energy += (fu-fv).norm2()*(*kuvMap)[e];				//Add edge contribution as kuv*(f(u)-f(v))^2
		else
			energy += (fu-fv).norm2();							//For tuette, contribution is just (f(u)-f(v))^2
	}
	return(energy);
}

//Computes laplacian for a vertex
Point calLaplace(Vertex* v,std::map<Edge*,double>* kuvMap){
	Point laplacian;
	HalfEdge *he = v->halfedge();
	do{											//Traverse neighbours
		Point fi = v->point();
		Point fj = he->source()->point();
		double kuv = (*kuvMap)[he->edge()];
		laplacian = laplacian + (fj-fi)*kuv;	//Add edge contibution as kuv*(f(u)-f(v))
		he = he->clw_rotate_about_target();
	}while(he != v->halfedge());
	return(laplacian);
}

//Computes updated map given relevant arguments using spherical constraints from formula used in Lai et al.
Point updatePt(Point pt, Point laplace, double step, double fh, double f2, double h2){
	double denom,num_a,num_b;									//we compute Y(t)=alpha(t)F + beta(t)H
	denom = (step*0.5);											//H=-laplacian
	num_a = (1.0 + (denom*fh));									//and t is step size
	denom = denom*denom;										//To simplify, we compute the numerators and denominator separately
	num_a = num_a*num_a - (denom*f2*h2);						//num_a = (1+(t/2)*(F.H))^2+(t/2)^2*F^2*H^2
	denom = 1.0 - (denom*fh*fh) + (denom*f2*h2);				//denom = 1-(t/2)^2*(F.H)^2+(t/2)^2*F^2*H^2
	num_b = -step*f2;											//num_b = -tF^2
	return(pt*(num_a/denom) + laplace*(num_b/denom)*(-1.0));	//Y(t) = (num_a/denom)F + (num_b/denom)H
}

int main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	std::ifstream in(argv[1]);
	if( argv[1][strlen( argv[1] )-1] == 'm'){
		MFileReader mf;
		mf.readToSolid(&mesh, in);
	}
	else{
		OBJFileReader of;
		of.readToSolid(&mesh, in);
	}
	mesh.UpdateNormals();											//Use library to compute normals

	/******************* Conformal Mapping Code *********************/

	std::map<Edge*,double> kuvMap;
	std::map<Vertex*,std::pair<Point*,Point*>> bbvarMap;			//Buffer structure to store data across iterations
																	//Needed for BB step size
	//Precompute Kuv for each edge
	SolidEdgeIterator eiter(&mesh);
	for(;!eiter.end();++eiter){
		Edge *e = *eiter;
		kuvMap[e] = calKuv(e);
	}

	//Initailize Conformal Mapping as Gauss Map
	SolidVertexIterator viter(&mesh);
	std::map<Vertex*,Point*> oldPointMap;
	//int vertexCount = 0;
	for(;!viter.end();++viter){
		Vertex *v = *viter;
		//++vertexCount;
		oldPointMap[v] = new Point();
		*(oldPointMap[v]) = v->point();
		v->point() = v->normal();									//Use vertex normals as initial map
	}

	double energy, deltaE;
	double step = 0.01;												//Set initial step size
	energy = calEnergy(&eiter,&kuvMap);
	deltaE = energy;
	std::cout<<"Initial Energy = "<<energy<<std::endl;

	//When using Barzilai–Borwein step size, required parameters will not exist in first iteration
	//So logic for first iteration differs from subsequent iterations
	//This is why it is outside the main update loop
	deltaE = energy;
	for(viter.reset();!viter.end();++viter){
		Vertex *v = *viter;
		Point laplacian;
		laplacian = calLaplace(v,&kuvMap);
		double fh,f2,h2;
		fh = (v->point())*((laplacian)*(-1.0));						//fh = <F(x),H(x)> where H(x)=-laplacian
		f2 = (v->point()).norm2();									//f2 = ||F(x)||^2
		h2 = laplacian.norm2();										//h2 = ||H(x)||^2
		Point newPt;
		newPt = updatePt(v->point(),laplacian,step,fh,f2,h2);
		std::pair<Point*,Point*> bbVars;
		bbVars.first = new Point();
		*(bbVars.first) = v->point();								//Save old F(x)
		bbVars.second = new Point();
		*(bbVars.second) = laplacian*f2 - v->point()*fh;			//Save old A(x)F(x)
		v->point() = newPt;											//F[k+1](x) = Y[k](x) ie Update map
		bbvarMap[v] = bbVars;
	}

	energy = calEnergy(&eiter,&kuvMap);
	std::cout<<"Harmonic Energy = "<<energy<<std::endl;
	deltaE -= energy;

	std::map<Vertex*,Point*> laplaceMap;							//Store laplacians in a buffer
	for(viter.reset();!viter.end();++viter){
		Vertex *v = *viter;
		laplaceMap[v] = new Point();
	}

	double epsilon=1e-10,rho=1e-4,delta=0.1,xi=0.85,ck=energy,qk=1;	//Hyperparameters for BB stepsize

	while(abs(deltaE) > epsilon){
		double agr_num=0.0,agr_denom=0.0;
		//This loop computes laplacians and aggregates for BB stepsize used in the next loop
		for(viter.reset();!viter.end();++viter){
			Vertex *v = *viter;
			*(laplaceMap[v]) = calLaplace(v,&kuvMap);
			Point dk_1,wk_1,afk;									//BB stepsize variables
			double f2,fh;
			fh = (v->point())*((*(laplaceMap[v]))*(-1.0));
			f2 = (v->point()).norm2();
			dk_1 = (v->point() - *(bbvarMap[v].first));				//D[k-1] = F[k]-F[k-1]
			agr_num += dk_1.norm2();								//Aggregate ||D[k-1]||^2
			afk = (*(laplaceMap[v]))*f2 - v->point()*fh;			//Get current A(x)F(x)
			wk_1 = afk - (*(bbvarMap[v].second));					//W[k-1] = A[k]F[k]-A[k-1]F[k-1]
			agr_denom += dk_1*wk_1;									//Aggregate D[k-1].W[k-1]
			*(bbvarMap[v].first) = v->point();						//Save currect F(x)
			*(bbvarMap[v].second) = afk;							//Save current A(x)F(x)
		}
		double stepk_1,val,new_energy;
		stepk_1 = agr_num/abs(agr_denom);							//Calculate t[k,1] from BB formula
		step = stepk_1;												//Set t[k] as t[k,1]
		do{															//Loop till BB condition is satisfied
			for(viter.reset();!viter.end();++viter){				//Calculate new map using current stepsize
				Vertex *v = *viter;
				double fh,f2,h2;
				fh = (*(bbvarMap[v].first))*((*(laplaceMap[v]))*(-1.0));
				f2 = (*(bbvarMap[v].first)).norm2();
				h2 = (*(laplaceMap[v])).norm2();
				v->point() = updatePt((*(bbvarMap[v].first)),(*(laplaceMap[v])),step,fh,f2,h2);
			}
			val = ck + rho*step*deltaE;								//Calculate condition parameters
			step *= delta;
			new_energy = calEnergy(&eiter,&kuvMap);
		}while(new_energy > val);									//Check E(Y[k](t[k])) > C[k] + rho*t[k]*E'(Y[k](0))
		ck = ((xi*qk*ck) + new_energy)/(xi*qk + 1);					//Compute C[k+1]
		qk = xi*qk + 1;												//Compute Q[k+1]
		std::cout<<"Harmonic Energy = "<<new_energy<<std::endl;
		deltaE = energy - new_energy;								//Compute change in energy
		energy = new_energy;
	}

	kuvMap.clear();
	
	std::ofstream os(argv[2]);
	
	if( argv[2][strlen( argv[2] )-1] == 'm'){
		// Write out the resultant .m file
		for(viter.reset();!viter.end();++viter){
			Vertex *v = *viter;
			os<<"Vertex "<<v->id();
			for (int i=0;i<3;++i)
				os<<" "<<(*(oldPointMap[v]))[i];
			os<<" {conformal=("<<(v->point())[0]<<" "<<(v->point())[1]<<" "<<(v->point())[2]<<")}"<<std::endl;
		}

		SolidFaceIterator fiter(&mesh);
		for(;!fiter.end();++fiter){
			Face *f = *fiter;
			os<<"Face "<<f->id();
			FaceVertexIterator fviter(f);
			for(; !fviter.end(); ++fviter)
			{
				Vertex *v = *fviter;
				os<<" "<<v->id();
			}
			os << std::endl;
		}
	}
	else{
		// Write out the resultant obj file
		int vObjID = 1;
		std::map<int, int> vidToObjID;

		std::ofstream os(argv[2]);

		SolidVertexIterator iter(&mesh);

		for(; !iter.end(); ++iter)
		{
			Vertex *v = *iter;
			Point p = v->point();
			os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			vidToObjID[v->id()] = vObjID++;
		}
		os << "# " << (unsigned int)mesh.numVertices() << " vertices" << std::endl;

		float u = 0.0, v = 0.0;
		for(iter.reset(); !iter.end(); ++iter)
		{
			Vertex *vv = *iter;
			std::string key( "uv" );
			std::string s = Trait::getTraitValue (vv->string(), key );
			if( s.length() > 0 )
			{
				sscanf( s.c_str (), "%f %f", &u, &v );
			}
			os << "vt " << u << " " << v << std::endl;
		}
		os << "# " << (unsigned int)mesh.numVertices() << " texture coordinates" << std::endl;

		SolidFaceIterator fiter(&mesh);
		for(; !fiter.end(); ++fiter)
		{
			Face *f = *fiter;
			FaceVertexIterator viter(f);
			os << "f " ;
			for(; !viter.end(); ++viter)
			{
				Vertex *v = *viter;
				os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
			}
			os << std::endl;
		}
	}

	os.close();
	return 0;
}