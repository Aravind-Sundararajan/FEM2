#include "PhyNode.h"
#include "PhyGlobal.h"

ostream& operator<<(ostream& out, const PhyNode& node)
{
	out << node.id;
	for (int i = 0; i < node.coordinate.size(); ++i)
	out  << '\t' << node.coordinate(i);
	// values
	out << '\n';
	for (int i = 0; i < node.nndof; ++i)
	out << node.ndof[i].v << '\t';

	// force
	out << '\n';
	for (int i = 0; i < node.nndof; ++i)
	out << node.ndof[i].f << '\t';
	//	out << '\n';

	if (verbose == true)
	{
		// position
		for (int i = 0; i < node.nndof; ++i)
		out << node.ndof[i].pos << '\t';
		out << '\n';

		// prescribed_boolean
		for (int i = 0; i < node.nndof; ++i)
		out << node.ndof[i].p << '\t';
	}
	return out;
}


void PhyNode::set_nndof(int nndofIn)
{
	nndof = nndofIn;
	ndof.resize(nndof);
}

void PhyNode::UpdateNodePrescribedDofForces(VECTOR& Fp)
{

	// complete

	//nodes
		for (int dofi = 0;dofi < nndof; ++dofi) //num dof for node (n)
		{
			if (ndof[dofi].p == true) //prescribed dof
			{
				posn = ndof[dofi].pos; //position of dof in global prescribed force F p
				ndof[dofi].f = Fp(-posn);
				//1. set prescribed dof force to corresponding force in global Fp (F p )
				//2. posn < 0; prescribed dof
			}
		}
	/*
	*/
}
