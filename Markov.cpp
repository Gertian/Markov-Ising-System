#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class IsingSystem;
ostream& operator<<(ostream& os,const IsingSystem& sys);

//this vector will hold the possible values for exp

//This class holds the ising system as an vector<vector<bool>> -> bool = true or false
class IsingSystem{
	////////////////
	//PRIVATE FIELDS
	////////////////
	private:
	vector<vector<int>> _system; // the system
	unsigned int _size;//the size of the system		

	double _J; //coupling from spin spin interactions/kT
	double _h; //external magnetic field*coupling/kT

	double _mag; //magnetisation of the system in units of 1 spin magnetisation
	double _energy; //energy of the system in units 1/kT
	
	bool _periodic;
	//////////
	//SETTERS
	/////////
	public:
	void setSpin(unsigned int row, unsigned int column, int a){
		this->_system[row][column] = a;
	}
	void addMag(double a){
		this->_mag += a;
	}	
	void addEnergy(double a){
		this->_energy += a;
	}

	///////////
	//GETTERS
	///////////
	public: 
	int getSpin(unsigned int row, unsigned int column)const{
		return this->_system[row][column];
	}
	unsigned int getSize()const{
		return this->_size;
	}
	double getJ()const{
		return this->_J;
	}
	double getH()const{
		return this->_h;
	}		
	double getMag()const{
		return this->_mag;
	}
	double getEnergy()const{
		return this->_energy;
	}

	//THis function will take one object in the system als argument (row, colomn) and return the amount of neighbours have the same spin
	//this is needed for the calculation of the potential energy of the system
	public:
	unsigned int getSameneighbours(unsigned int row, unsigned int column)const{
		//init the returnvalue at -1 in order to correct for the element-element count that will be made...
		unsigned int returnvalue = -1;
		//first we retrieve the spin of the current object
		int centralspin = getSpin(row, column);
		for(int i = -1; i<=1; i++ ){
			for(int j = -1; j<=1; j++){
				int thisspin;
				if(!_periodic){
					//we will need to check that we are not out of bounds !
					if(row + i >= 0 && row + i < _size && column + j >= 0 && column + j < _size){
						//in this case we are safely in the system and will get no OOB error
						thisspin = getSpin(row + i, column + j);
					}
				}else{
					//now we will take it into account even when we are OOB. to always get the copy inside the system
					//we simply do row+i+_size%_size idd for column.
					thisspin = getSpin((row + i + _size)%_size, (column+j+_size)%_size);
				}
				if(thisspin == centralspin){
					returnvalue++;
				}
			}
		}
		return returnvalue;
	}	
	
	//THis function retursn the total amount of neighbours for that particle. 
	public:
	unsigned int getAmountneighbours(unsigned int row, unsigned int column)const{
		if(_periodic){
			return 8;
		}
		//check if the particle is on the left side or right side
		unsigned int neighbours = 8;
		bool a = false;
		bool b = false;
		if(row == 0 || row == _size-1){
			a = true;
		}
		if(column == 0 || column == _size-1){
			b = true;
		}
		if(a == true || b == true){
			neighbours = 5;
			if(a == true && b == true){
				neighbours = 3;
			}
		}
		return neighbours;
	}

	//this function calculates the total spin in the system
	//this is needed for the init of the system
	public:
	int calcSpin()const{
		int spin = 0;
			for(unsigned int i = 0;i<_size; i++){
				for(unsigned int j =0; j<_size; j++){
					spin += getSpin(i,j);
				}
		}
		return spin;
	} 


	//this function calculates the total energy in the system
	//this is needed for the init of the system
	public:
	double calcEnergy()const{
		double energy = 0;
		//loop over all the particles in the system
		for(unsigned int i = 0; i< _size; i++){
			for(unsigned int j = 0; j<_size; j++){
				//add the same spin contribution of the particle
				energy -= 0.5*_J*getSameneighbours(i,j);
				//add the contribution of the spin/field interation of that particle
				energy -= _h*getSpin(i,j);
			}
		}
		return energy;
	}

	//this function calculates the energy change in the system if i,j were to be flipped
	public:
	double calcEnergydifference(unsigned int row, unsigned int column)const{
		//energy b4 = -_J*(-{neigbours - sameneighbours} + sameneighbours ) + spincont
		//          = -_J*(+2*SN - N)
		//energy a  = -_J*(+{neigbours - sameneighbours} - sameneighbours ) - spincont
		//          = -_J*(-2*SN + N)	
		// a - b4   = -_J*(-4*SN + 2*N)

		//TODO This formula is not reversible ! It should be adapted ?! 
		double a = -this->_J*(-4.*this->getSameneighbours(row, column) + 2.*this->getAmountneighbours(row, column)) + 2*this->getSpin(row, column)*this->_h;
		return a;
	}

	public:
	//constructor of the Isingsystem class. The arguments are:
	//->a systemsize, amount of particles will be this squared!
	//->a coupling constant(in units of inverse temperature*kB). How hard is it to flip one spin ?!
	//->All particles will be initialised on the true value!
	IsingSystem(unsigned int  size, double J,double h, bool periodic){
		for(unsigned int row = 0; row < size; row++){
			vector<int> a(size, 1);
			_system.push_back(a);
		}
		_size = size;//the size of the system		
	
		_J = J; //coupling from spin spin interactions/kT
		_h = h; //external magnetic field*coupling/kT
	
		_mag = calcSpin(); //magnetisation of the system in units of 1 spin magnetisation
		_energy = calcEnergy(); //energy of the system in units 1/kT

		_periodic = periodic;
	}

	//This function will get an i,j location. Get the amount of same spin friends. 
	//calculate the odds of that spin changing using the metropolis principle
	//generate a random number and see if it actually spins
	//returns true if spin, false if not spint
	public:
	bool getsSpinned(unsigned int row, unsigned int column)const{
		//first we calculate the exp(-Benergy difference)
		double diff = calcEnergydifference(row, column);
		double probability  = exp(-diff);
		//generate a random number from zero to one
		double random = (double)rand()/RAND_MAX;
		if(probability >= random ){
			return true;
		}else{
			return false;
		}
	}

	//this function will get an i,j location and flip the spin
	void flipSpin(unsigned int row, unsigned int column){
		//adjust the energy
		this->addEnergy(this->calcEnergydifference(row, column));
		//adjust the magnetisation
		this->addMag(-2*this->getSpin(row, column));
		//adjust the actual spin
		this->setSpin(row, column, -this->getSpin(row, column));
	}

	//this function will make one markov step --> check if one has to be spinned and spin it accordingly
	void MarkovStep(){
		//loop over the system
		for(unsigned int i =0; i<this->getSize(); i++){
			for(unsigned int j = 0; j<this->getSize(); j++){
				if(this->getsSpinned(i, j)){
					this->flipSpin(i,j);
				}
			}
		}
	}
	
	//this fucntion wil excectute n markov steps for a given system
	void MarkovSteps(unsigned int amount){
		ofstream spins("Results.md");
		ofstream data("energyandmag.md");
		cout << "Simulating for " << amount << " MC steps" << endl;
		cout << "-------------------------------------------------------------------" << endl;
		spins << *this << endl << endl;
		data  << 0 << "\t" << this->getEnergy() << "\t" << this->getMag()/this->getSize()/this->getSize() << endl;
		for(unsigned int i = 1; i <= amount; i++){
			this->MarkovStep();
			cout << "Progres: " << i << '\r';
			cout.flush();
	
			spins << *this << endl << endl;
			data  << i << "\t" << this->getEnergy() << "\t" << this->getMag()/this->getSize()/this->getSize() << endl;
		}
		cout << "Simulation completed" << endl;
		cout << "Making the plots" << endl;
		cout << "-------------------------------------------------------------------" << endl;
		string a = "gnuplot -e \"amount = "+to_string(amount)+"\" Plot.gnu";
		system(a.c_str());
		cout << endl;
		cout << "Making a call to ffmpeg to make a movie out of the simulation" << endl;
                cout << "-------------------------------------------------------------------" << endl;
                system("ffmpeg -y -framerate 5 -i Plot%03d.png -c:v libx264 -r 30 -pix_fmt bgr565 out.mp4");
                system("rm Plot*.png");
                cout << endl << endl <<"done" << endl;
	}


};

ostream& operator<<(ostream& os,const IsingSystem& sys){
	double size = sys.getSize();
	for(unsigned int i = 0; i<size; i++){
		for(unsigned int j = 0; j<size; j++){
			os << i/size << "\t" << j/size << "\t" << sys.getSpin(i,j) << endl;
		}
	}
	return os;
}

int main(){
	unsigned int size;
	double J;
	double h;
	unsigned int amount;	
	bool periodic;

	cout << "What should the size of the system be ? size >= 3 for valid results" << endl;
	cin >> size;	
	if(size > 5000){
		cout << "THis system size will likely crash the system" << endl;
		cout << "Killing progres for safety" << endl;
		exit(0);	
	}

	cout << "What should the spin-spin coupling be ? This is in units of coupling*single spin field/kT, >0" << endl;
	cin >> J;	

	cout << "What should the spin-field coupling be ? This is in units of coupling*external field/kT" << endl;
	cin >> h;

	cout << "how many Mc steps do you want to do ?" << endl;
	cin >> amount;

	cout << "Do you want to impose periodic boundary conditions? " << endl;
	cin >> periodic;

	IsingSystem sys(size, J, h, periodic);
		
	sys.MarkovSteps(amount);
}





