#include "CFEMTypes_Global.h"
#include "FEMSolver.h"
#include "PhyElement.h"
#include "LAFuncs.h"
#include "PhyGlobal.h"
#include "CFEMTypes_Global.h"


// in C++ do not write friend again (similar to virtual)
void FEMSolver::Input(istream& in)
{
	// READING Nodes, ... OMITTED
	// ....
	// ....
	// ....
	string buf;
	in >> buf >> dim; // dim 2
	in >> buf >> ndofpn; // ndofpn 2
	in >> buf >> buf >> nNodes; // Nodes nNodes 3
	in >> buf >> buf; //id crd
	nodes.resize(nNodes);
	for (int i = 0; i <nNodes; ++i)
	nodes[i].set_nndof(ndofpn);

	int tmpi;
	for (int i = 0; i < nNodes; ++i)
	{
		in >> tmpi;
		if (tmpi != (i + 1))
		{
			THROW("incorrect id")
		}
		// use tmpi - 1 = i for id instead to simplipy accessing nodes
		// for a more robust implementation that you can have any node ids, read with arbitrary
		nodes[i].id = i;
		nodes[i].coordinate.resize(dim);
		for (int j = 0; j < dim; ++j)
		in >> nodes[i].coordinate(j);
	}
	in >> buf >> buf >> ne;
	in >> buf >> buf >> buf >> buf >> buf;
	pes.resize(ne);
	ElementType eType;
	int matID;
	int nNodeInElement;
	PhyElement* pe;

	for (int i = 0; i < ne; ++i)
	{
		in >> tmpi;
		if (tmpi != (i + 1))
		THROW("incorrect id");

		in >> tmpi;
		eType = (ElementType)tmpi;

		pes[i] = PhyElementFactory(eType);
		pes[i]->id = tmpi;

		// 1. we don't need to use pe instead of pes[i]. It just makes read and write simpler
		// 2. we can use pe instead pes[i] because it's a pointer (it's an address).
		// 3. OBVIOUSLY we cannot do this trick with nonpointer data

		pe = pes[i];
		in >> pe->matID;		// longer way which was fine: in >> pes[i]->matID;
		// another way not recommended (*pe).matID
		// ptr			*ptr	object
		// object		&object	address of the object

		//int nNodeInElement;
		in >> nNodeInElement;
		vector<int> eNodesTmp(nNodeInElement);
		vector <PhyNode*> eNodePtrsTmp(nNodeInElement);
		for (int j = 0; j < nNodeInElement; ++j)
		{
			in >> eNodesTmp[j];
			--eNodesTmp[j];
			eNodePtrsTmp[j] = &nodes[eNodesTmp[j]]; // safe here because nodes size never is going to change. If not this causes a very nasty bug to fix ...
		}
		pe->setNodeConnectivity_Sizes(nNodeInElement, ndofpn, eNodesTmp, eNodePtrsTmp);
	}
	in >> buf >> buf >> np;
	in >> buf >> buf >> buf; //node  node_dof_index  value

	int nodeid;
	int dofid;
	double value;
	PhyDof* dofPtr;
	for (int i = 0; i < np; ++i)
	{
		in >> nodeid >> dofid >> value;
		--nodeid;
		--dofid;
		// could have daone the last three as the following -- and ++ after the parameter does
		//	in >> nodeid-- >> dofid-- >> value;

		// good practice to with a shorter pointer rather than the full name
		dofPtr = &nodes[nodeid].ndof[dofid]; // & to get the point
		dofPtr->p = true;
		dofPtr->v = value;
		// need to
		// A. assign these values to nodes (Step 4 in course notes)
	}

	int nnzdof; // num of nonzero force free Dofs;
	in >> buf >> buf >> nnzdof;

	in >> buf >> buf >> buf ; //node  node_dof_index  value
	for (int i = 0; i < nnzdof; ++i)
	{
		in  >> nodeid >> dofid >> value;
		--nodeid;
		--dofid;

		// good practice to with a shorter pointer rather than the full name
		dofPtr = &nodes[nodeid].ndof[dofid]; // & to get the point
		//		dofPtr->p = false; // no need for this (default is false)
		dofPtr->f = value; // force is given
	}

	in >> buf >> buf >> nmats;
	in >> buf >> buf >> buf; // id  numPara  Paras

	int numParas, matid;
	for (int i = 0; i < nmats; ++i)
	{
		in >> matid >> numParas;
		//		 	--matid;
		// 			if (matid !=i)
		//				THROW("wrong material id\n")

		mats[matid].setSize(numParas);
		for (int j = 0; j < numParas; ++j)
		in >> mats[matid].paras(j);
	}

	for (int e = 0; e < ne; ++e)
	{
		pe = pes[e];
		pe->setGeometry();
		matID = pe->matID;
		pe->setInternalMaterialProperties(&mats[matID]);
	}
	//return in;
	//verifying the input
	cout << "dim " << dim << "\n";
	cout << "ndofpn " << ndofpn << "\n";
	cout << "nNodes " << nNodes << "\n";
	cout << "nElements " << ne << "\n";
}

istream& operator>>(istream& input, FEMSolver& dat)
{
	dat.Input(input);
	return input;
}

ostream& operator<<(ostream& out, const FEMSolver& dat)
{
	out << "Nodes\n";
	out << "nNodes\t" << dat.nNodes << '\n';
	out << "id\tcrd\n";
	out << "values\nforces\n";
	if (verbose == true)
	{
		out << "position(verbose)\n";
		out << "prescribed_boolean(verbose)\n";
	}
	for (int node = 0; node < dat.nNodes; ++node)
	out << dat.nodes[node] << '\n';
	// Complete the function
	out << "Elements\n";
	out << "ne\t" << dat.ne << "\n";
	out << "id ElementType\n";
	out << "forces(verbose)\n";
	out << "specific output\n";

	for (int e=0; e<dat.ne; ++e)
	out << (*dat.pes[e]) << '\n';

	return out;
}

FEMSolver::FEMSolver(int dimIn)
{
	dim = dimIn;
}

FEMSolver::~FEMSolver()
{
	// ne should be equal to pes.size()
	// still a better practice is
	for (int i = 0; i < pes.size(); ++i)
	delete pes[i];
}


void FEMSolver::FEMSolve(string& runName, bool verboseIn)
{
	verbose = verboseIn;
	string inputFileName;
	inputFileName = runName + ".txt";
	fstream in(inputFileName.c_str(), ios::in);
	if (in.is_open() == false)
	{
		cout << "input file name\t" << inputFileName << " does not exist\n";
		THROW("file does not exist\n");
	}
	// reading data
	Input(in);
	// can do it as
	//	in >> (*this);
	in.close();



	/////////////////////////////////////////////////////////////////////////
	// steps

	// Step 3
	setSizes();
	cout << "1 \n";
	// Step 4; set prescribed dofs: already done when reading the input file
	// Step 5: Set global free nodal dof: already done when reading the input file
	// Step 6 and Step 7: dof positions; Step 7: Set F
	setPositions_F();
	cout << "2 \n";
	// Step 8: Element dof maps Me
	// Step 9: Set element dofs ae
	setElementDofMap_ae();
	cout << "3 \n";
	// Step 10: Compute element stiness
	//Calculate_ElementStiffness_Force();

	// Step 11: Assembly from local to global system
	Assemble();
	cout << "4 \n";
	// Step 12: Solve global (free) dof a from Ka = F
	// successful solution returns true
	if (Solve_Dofs() == false)
	THROW("Matrix solve failed\n");
	// Step 13: Assign a to nodes and elements
	Assign_dof();
	cout << "5 \n";
	// Step 14: Compute prescribed dof forces
	UpdateFpNodalPrescribedForces();
	cout << "6 \n";
	/////////////////////////////////////////////////////////////////////////
	// output
	string outputFileName;
	outputFileName = runName + "Output.txt";
	fstream out(outputFileName.c_str(), ios::out);
	out << (*this);

	cout << "wrote to file \n";
}

void FEMSolver::setSizes()
{

	ndof = nNodes*ndofpn;
	nf = ndof - np;
	// Complete
	K.resize(nf,nf);
	F.resize(nf);
	Fp.resize(nf);
	//PhyElement* pe;
	//for (int e = 0; e < ne; ++e)
	//pes[e]->setNodeConnectivity_Sizes(nNodeInElement, ndofpn, eNodesTmp, eNodePtrsTmp);
}

void FEMSolver::setPositions_F()
{
	// Complete
	int posf = 0;
	int posp = 0;
	for (int n = 0;n < nNodes; ++n) {
		for (int dofi = 0; dofi <ndofpn; ++dofi){ // num dof for node (n)
			cout << "d is: " <<dofi << "\n";
			if (nodes[n].ndof[dofi].p == true){ // prescribed dof
				posp = posp - 1;
				nodes[n].ndof[dofi].pos = posp;
			}
			else {  //free dof
				posf = posf + 1;
				nodes[n].ndof[dofi].pos = posf;
				F(posf-1) = nodes[n].ndof[dofi].f;
			}
		}
	}
cout << "the F positions is :  "<< F<< "\n";
}


void FEMSolver::setElementDofMap_ae()
{
	cout << "trying to set elementDofMap \n";
	PhyElement* pe;
	for (int e = 0; e < ne; ++e)
	{
		cout << "We are looping over element" << e << "\n";
		pes[e]->setElementDofMap_ae(ndofpn);
	}
}


void FEMSolver::Calculate_ElementStiffness_Force()
{
	PhyElement* pe;
	for (int e = 0; e < ne; ++e)
	pes[e]->Calculate_ElementStiffness_Force();
}

void FEMSolver::Assemble()
{
	PhyElement* pe;

	cout << "Assembly: ne is " << ne << "\n";
	for (int e = 0; e < ne; ++e)
	{
	cout << "e is" << e << "\n";
	pes[e]->AssembleStiffnessForce(K, F);
	}
}

bool FEMSolver::Solve_Dofs()
{
	dofs.resize(nf);
	dofs = F;
	cout << "K\n" << K << endl;
	int isNonsingular;
if (verbose)
{
	// to save K for verbose Output
	MATRIX Kbk;
	Kbk = K;
	isNonsingular = !LUsolve(K, dofs);
}
else
{
isNonsingular = !LUsolve(K, dofs);
}
return (isNonsingular != 0);
}

void FEMSolver::Assign_dof()
{
	int posn = 0;
	//Complete
	//assign to nodes
	for (int n = 0; n < nNodes; ++n)
	{
		for (int dofi = 0; dofi< nodes[n].nndof; ++dofi)//num dof for node (n)
		{
			if (nodes[n].ndof[dofi].p == false) //free dof
			{
				posn = nodes[n].ndof[dofi].pos; //position of dof in global free F
				nodes[n].ndof[dofi].v = dofs(posn-1); //set free dof val to corresponding val in global dofs (a)
			}
		}
	}
	// assign to elements
	for (int e = 0; e < ne; ++e){// loop over elements
		for (int i =0; i < pes[e]->nedof; ++i){ //loop over element dofs; nedof = # dof (nedof )
			posn = pes[e]->dofMap[i]; //corresponding global position using dofMat (Met)
			if (posn > 0) //free dof
			{
				pes[e]->edofs(i) = dofs(posn-1);
				//set free element dof ae to corresponding val in global dofs (a)
			}
		}
	}

}

void FEMSolver::UpdateFpNodalPrescribedForces()
{
	// Complete
	PhyElement* pe;
	for (int e = 0; e < ne; ++e) //loop over the elements
	pes[e]->UpdateElementForces_GlobalFp(F);

	for (int n = 0; n < nNodes; ++n) //loop over the elements
	nodes[n].UpdateNodePrescribedDofForces(F);
}
/*
*/
