#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* definir a fonction template pour calculer le produit interne
   entre deux valarray
   inputs:
     array1: (valarray<T>)(N) vecteur de taille N
     array2: (valarray<T>)(N) vecteur de taille N
   outputs:
     produitInterne: (T) produit entre les deux vecteurs
*/


template<typename T> T scalarProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
} 

/* definir a fonction template pour calculer la norm2 d'un valarray
   inputs:
     array: (valarray<T>)(N) vecteur de taille N
   outputs:
     norm2: (T) norm2 du vecteur
*/
template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
} 

/* definir a fonction template pour calculer le produit vecteur
   entre 2 valarray de dimension 3
   inputs:
     array1, array2: (valarray<T>)(N) vecteurs de taille N
   outputs:
     produitVecteur: (T) produit vectoriel array1 x aray2 
*/
template<typename T> T vectorProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // initialiser le nouveau valarray
  valarray<T> array3=valarray<T>(3);
  // calculer le produit vecteur
  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premier composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante
  // retourner le array3
  return array3;
} 


/* La class Engine est le moteur principal de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{

protected:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  // definition des variables
  unsigned int nsteps=100000; // Nombre de pas de temps
  double rmax = 120.0e9;
  double G= 6.67430e-11;
  double m1=1.989e30;      // mass du soleil (ou masse du corps central)
  double m2=685;
  double m3=10e-31;
  double tfin=9400753.64319;      // Temps final
  double epsilon=10000;
  
  
  valarray<double> Y=valarray<double>(0.e0,12); // position through time 
  valarray<double> Y0=valarray<double>(0.e0,12); // initial position
  
  
  
  unsigned int sampling=1; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  // TODO calculer l'energie mecanique et le moment magnetique
  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
       write: (bool) ecriture de tous les sampling si faux
 
  */  
  
  double energy(){
	  double energy;
	  energy = (1.0/2.0)*m2*(Y[8]*Y[8]+Y[9]*Y[9])-G*m1*m2/(sqrt( pow(Y[2]-Y[0],2) + pow(Y[3]-Y[1],2) ));
	  return energy;
  }
  
  
  void printOut(bool write)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << t << " " << Y[0] << " " << Y[1] << " " << Y[2] << " "<< Y[3] << " " << Y[4] << " " << Y[5] << " " 
      << Y[6] << " " << Y[7] << " " << Y[8] << " " << Y[9] << " " << Y[10] << " " << Y[11] << " " << energy() << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;

//---------------- protected: ------------------------------------------

  // donnes internes
  double t,dt;        // Temps courant pas de temps

   /* Cette methode calcule le champ magnetique en y
     outputs:
       By: (double) champ magnetique en y
  */ 

  // TODO
  /* Cette methode calcule l'acceleration
     outputs:
       a: (valarray<double>)(2) acceleration dans les directions (x,z)
  */
  valarray<double> f(valarray<double> y) const // l'accéleration est déduite de la 2eme loi de Newton comme décrit dans le rapport
  {
    
    //normes 
    // norme de r1-r2 puis 3
    double r12 = pow(sqrt(( pow(y[0]-y[2],2) + pow(y[1]-y[3],2) )),3.0);
    
    // norme de r2-r3 puis 3
    double r13 = pow(sqrt(( pow(y[0]-y[4],2) + pow(y[1]-y[5],2) )),3.0);
    
    // norme de r1-r3 puis 3
    double r23 = pow(sqrt(( pow(y[2]-y[4],2) + pow(y[3]-y[5],2) )),3.0);
    
    valarray<double> a=valarray<double>(0.e0,12);		
    
    // compute the acceleration
    a[0]= y[6];
    a[1]= y[7];
    a[2]= y[8];
    a[3]= y[9];
    a[4]= y[10];
    a[5]= y[11];
    
    a[6]= -G*( m2*(y[0]-y[2])/(r12) +m3*(y[0]-y[4])/(r13) );
    a[7]= -G*( m2*(y[1]-y[3])/(r12) +m3*(y[1]-y[5])/(r13) );
    
    a[8]= -G*( m1*(y[2]-y[0])/(r12) +m3*(y[2]-y[4])/(r23) );
    a[9]= -G*( m1*(y[3]-y[1])/(r12) +m3*(y[3]-y[5])/(r23) );
    
    a[10]= -G*( m1*(y[4]-y[0])/(r13) +m2*(y[4]-y[2])/(r23) );
    a[11]= -G*( m1*(y[5]-y[1])/(r13) +m2*(y[5]-y[3])/(r23) );
    

    
    
    
    
    return a;
  }

public:

  /* Constructeur de la classe Engine
     inputs:
       configFile: (ConfigFile) handler du fichier d'input
  */
  Engine(ConfigFile configFile)
  {
    // variable locale
    
    // Stockage des parametres de simulation dans les attributs de la classe
    /*tfin     = configFile.get<double>("tfin");		 // lire la temps totale de simulation
    nsteps   = configFile.get<unsigned int>("nsteps");   // lire la nombre de pas de temps
    E0       = configFile.get<double>("E0");		 // lire l intensite du champ electrique
    x0[0]    = configFile.get<double>("x0");		 // lire composante x position initiale             à remplacer !!
    x0[1]    = configFile.get<double>("z0");		 // lire composante z position initiale
    v0[0]    = configFile.get<double>("vx0");		 // lire composante x vitesse initiale
    v0[1]    = configFile.get<double>("vz0");		 // lire composante z vitesse initiale
    sampling = configFile.get<unsigned int>("sampling"); // lire le parametre de sampling*/


    Y0[0] = configFile.get<double>("x1");
    Y0[1] = configFile.get<double>("y1");  
    Y0[2] = configFile.get<double>("x2");
    Y0[3] = configFile.get<double>("y2");
    Y0[4] = configFile.get<double>("x3");
    Y0[5] = configFile.get<double>("y3");
    Y0[6] = configFile.get<double>("v1x");
    Y0[7] = configFile.get<double>("v1y");
    Y0[8] = configFile.get<double>("v2x");
    Y0[9] = configFile.get<double>("v2y");
    Y0[10] = configFile.get<double>("v3x");
    Y0[11] = configFile.get<double>("v3y");
    
    m1 = configFile.get<double>("mass1");
    m2 = configFile.get<double>("mass2");
    m3 = configFile.get<double>("mass3");
    tfin     = configFile.get<double>("tfin");
    nsteps   = configFile.get<double>("nsteps");
    sampling = configFile.get<double>("sampling");
    epsilon = configFile.get<double>("tol");
    

    dt= tfin / nsteps;          // calculer le time step


    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str()); 
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  void run()
  {
    t = 0.e0; // initialiser le temps
    Y = Y0;   // initialiser le vecteur Y
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps
    /*while(abs(tfin-t)>dt/2.0) // boucle sur temps
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }*/
    
    
    for(unsigned int i(0); i<nsteps+2	; ++i) // boucle sur tout pas de temps
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  }

};




// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineRungeKutta4: public Engine
{
public:

  // construire la class Engine
  EngineRungeKutta4(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du 7ouvement en utilisant
     le scheme: Runge-kutta 2
  */
  void step()
  {

      // initialisations: 
      valarray<double> y=Y;
      
      valarray<double> k1=valarray<double>(0.e0,12);
      valarray<double> k2=valarray<double>(0.e0,12);
      valarray<double> k3=valarray<double>(0.e0,12);
      valarray<double> k4=valarray<double>(0.e0,12);
      
      
      
      
      
      k1 = dt * f(y);
      k2 = dt * f(y+(1.0/2.0)*k1);
      k3 = dt * f(y+(1.0/2.0)*k2);
      k4 = dt * f(y+k3);
      
      
      
      
      y = y+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
      t += dt; 
      
      Y=y; 
      
  } 
  
 
};



// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineTempsAdaptatif: public Engine {



public:

  // construire la class Engine
  EngineTempsAdaptatif(ConfigFile configFile): Engine(configFile) {}
 
  void step(){
	  double deltat = dt;
    // calcul du pas entier : 
      valarray<double> y1=Y;
      
      // initialisation
	 cout << ">>>" << dt << endl ;
      
      valarray<double> k1=valarray<double>(0.e0,12);
      valarray<double> k2=valarray<double>(0.e0,12);
      valarray<double> k3=valarray<double>(0.e0,12);
      valarray<double> k4=valarray<double>(0.e0,12);
      
      k1 = deltat * f(y1);
      k2 = deltat * f(y1+(1.0/2.0)*k1);
      k3 = deltat * f(y1+(1.0/2.0)*k2);
      k4 = deltat * f(y1+k3);
      y1 = y1+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
      

  
    // calcul du pas divisé en 2 : 
      valarray<double> y2=Y; 
      valarray<double> k21=valarray<double>(0.e0,12);
      valarray<double> k22=valarray<double>(0.e0,12);
      valarray<double> k23=valarray<double>(0.e0,12);
      valarray<double> k24=valarray<double>(0.e0,12);    
       
      k21 = (deltat/2.0) * f(y2);
      k22 = (deltat/2.0) * f(y2+(1.0/2.0)*k21);
      k23 = (deltat/2.0) * f(y2+(1.0/2.0)*k22);
      k24 = (deltat/2.0) * f(y2+k23);
      y2  = y2+(1.0/6.0)*(k21+2.0*k22+2.0*k23+k24);
      
      
      valarray<double> y3=y2; 
      valarray<double> k31=valarray<double>(0.e0,12);
      valarray<double> k32=valarray<double>(0.e0,12);
      valarray<double> k33=valarray<double>(0.e0,12);
      valarray<double> k34=valarray<double>(0.e0,12);     
      
      k31 = (deltat/2.0) * f(y3);
      k32 = (deltat/2.0) * f(y3+(1.0/2.0)*k31);
      k33 = (deltat/2.0) * f(y3+(1.0/2.0)*k32);
      k34 = (deltat/2.0) * f(y3+k33);
      y3  = y3+(1.0/6.0)*(k31+2.0*k32+2.0*k33+k34); 
    
    
    // comparaison des deux solutions : 
    valarray<double> deltay=y1-y3;
    double d = norm2(deltay);
    
    //+++++++++++++++ ++++++++++++++++
    if(d<epsilon) { 
		cout << d << endl;
		
		cout << ""<< endl;
		t = t+dt;
		dt = deltat*pow((epsilon/ (d+pow(10.0,-10) ) ),1.0/5.0);
	}
    else{
    do{ 
	    
	    //cout << dt << endl;
        cout << "avant " << d<< "," << deltat << endl ;
		deltat = 0.95*deltat*pow((epsilon/(d+pow(10.0,-10))),1.0/5.0); 
		
		y1=Y;
		k1 = deltat * f(y1);
		k2 = deltat * f(y1+(1.0/2.0)*k1);
		k3 = deltat * f(y1+(1.0/2.0)*k2);
		k4 = deltat * f(y1+k3);
		y1 = y1 + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
		
		
		y2=Y;
		k21 = (deltat*0.5) * f(y2);
		k22 = (deltat*0.5) * f(y2+(1.0/2.0)*k21);
		k23 = (deltat*0.5) * f(y2+(1.0/2.0)*k22);
		k24 = (deltat*0.5) * f(y2+k23);
		y2 = y2 + (1.0/6.0)*(k21+2.0*k22+2.0*k23+k24);
     
		y3=y2;
		k31 = (deltat*0.5) * f(y3);
		k32 = (deltat*0.5) * f(y3+(1.0/2.0)*k31);
		k33 = (deltat*0.5) * f(y3+(1.0/2.0)*k32);
		k34 = (deltat*0.5) * f(y3+k33);
		y3 =  y3 + (1.0/6.0)*(k31+2.0*k32+2.0*k33+k34); 
		
    
		
		deltay=y3-y1;
        d = norm2(deltay);
        cout << "après " << d << "," << deltat << endl ;
        cout << "while" << endl << endl;
    
	} while(d>epsilon);
	   dt = deltat;
	   t = t+dt;
    }
    

    cout << "----------------------"<< dt <<"----------------------" << endl;
    Y = y1;    
   }
	
};



// programme
int main(int argc, char* argv[])
{
  string inputPath("config.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique ("Euler"/"E", "EulerCromer"/"EC" ou "RungeKutta2"/"RK2")
  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser
  if(schema == "RungeKutta4" || schema == "RK4")
  {
    // initialiser une simulation avec schema runge-kutta 2
    engine = new EngineRungeKutta4(configFile);
  }
  else if(schema == "TempsAdaptatif" || schema == "TA")
  {
    // initialiser une simulation avec schema Boris Buneman
    engine = new EngineTempsAdaptatif(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}
