#include <iostream>
#include <vector>
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
 
typedef std::vector<boost::multiprecision::cpp_rational> stdvectorboostmultiprecisioncpp_rational;


int gcd(int a , int b){
  if (b==0)
  {
    return a;

  }
  if (a == 0)
  {
   return 1;
  }
  
  return gcd(b,a%b);
}

class rational {
  private:
    void simplify(rational &r){
      int  g;
      g = gcd(r.denominator,r.numerator);
      r.numerator = r.numerator/g;
      r.denominator = r.denominator/g;
      if (r.denominator <=0)
      {
        r.numerator = r.numerator *(-1);
        r.denominator = r.denominator *(-1);
      }
      
    }

  
   

  public:
    int numerator;
    int denominator;
    rational(){
      numerator = 0;
      denominator = 1;
    }
    rational(int n,int d){
      if (n>=0 && d>=0)
      {
        numerator = n;
        denominator = d;
      }
      if (n<=0 && d>=0)
      {
        numerator = n;
        denominator = d;
      }
      if (d<=0)
      {
        numerator = (-1)*n;
        denominator = (-1)*d;
      }
      
    }
  

 void show(){
    cout << this->numerator << "/" << this->denominator  ;
}
 rational operator+ (  const rational& rational_b){
   rational tmp;
   tmp.numerator  = (rational_b.numerator * this->denominator ) +(this->numerator * rational_b.denominator);
   tmp.denominator  = rational_b.denominator * this->denominator;
   simplify(tmp);
   
   return tmp;
 }
 rational operator- (  const rational& rational_b){
   rational tmp;

   tmp.numerator  = (this->numerator * rational_b.denominator) - (this->denominator *rational_b.numerator); 
   tmp.denominator  = rational_b.denominator * this->denominator;
   simplify(tmp);
   if (tmp.denominator <=0)
   {
    tmp.numerator  = (-1)*tmp.numerator;
    tmp.denominator  = (-1)*tmp.denominator;

   }
   
   return tmp;
 }
 rational operator*( const rational& rational_b){
  rational tmp;
  tmp.numerator  = this->numerator *rational_b.numerator;
  tmp.denominator = this->denominator * rational_b.denominator;
  simplify(tmp);

  return tmp;
 }
 
 rational operator/( const rational& rational_b){
  rational tmp;
  tmp.numerator  = this->numerator *rational_b.denominator;
  tmp.denominator = this->denominator * rational_b.numerator;
  
  simplify(tmp);

  return tmp;
 }

 
 rational& operator= (  const rational& rational_b){
    
    this->numerator = rational_b.numerator;
    this->denominator = rational_b.denominator;
  return *this;
}

rational& operator+=( const rational& rational_b){
   *this = *this + rational_b;

   return *this;
 }
};
class poly{
  private: 
    
    

  public:  
  rational coeff[10]; 
  long int deg ;  
   poly(){
    deg = 0;
    for (long int i = 0; i < 10 ; i++)
    { //rational b = {0,0};
      coeff[i].numerator = 0;
      coeff[i].denominator = 1;
    }
    
   }

  void set(rational a, long int d){
    coeff[d].numerator = a.numerator;
    coeff[d].denominator = a.denominator;
    deg = degree();
  }

  void set_to_zeros(){
    for (long int i = 0; i <= deg; i++)
    {
    coeff[i].numerator = 0;
    coeff[i].denominator = 1;
    }
    deg = degree();
  }

  long int degree()  {
    long int d = 0;
    for (long int  i = 0; i < 10; i++)
    {
      if (coeff[i].numerator != 0  ){ 
        d = i;
         
      }
    }
    deg = d;
    return d;
  }
 // poly poly_modular(poly poly_a , poly poly_b);
  poly operator= (const poly &poly_b){
      for (long int  i = 0; i <= max(poly_b.deg,this->deg); i++)
      {
        this->coeff[i].numerator = poly_b.coeff[i].numerator;
        this->coeff[i].denominator = poly_b.coeff[i].denominator;
      }
      deg = this->degree();
      return *this;
  }
  rational operator[](long int idx){
    return coeff[idx];
  }
};
typedef std::vector<rational> stdvectorrational;

struct FiveTuple: poly{
   // using sstdvectorrational::stdvectorrational;
   // typedef typename FiveTuple::value_type Entry;
    template <typename T>
        FiveTuple(T in) : poly() {
            
        }

    void zero(){
     for (long int i = 0; i <= deg; i++)
     {
     coeff[i].numerator = 0;
     coeff[i].denominator = 1;
     }
     deg = degree();
     }
    void one(){
      rational one(1,1);
      Zero();
      set(one, 0);
    }

    static FiveTuple Zero(){
      return FiveTuple(2);

    }
    static FiveTuple One(){
      rational one(1,1);
      FiveTuple p = FiveTuple(0);
      p.set(one,0);
      return p;

    }
    static FiveTuple sqrt2(){
      rational half(1,2);
      FiveTuple p = FiveTuple(0);
      p.set(half,1);
      p.set(half,3);

      return p;
    }
    static FiveTuple Rand();
    
    
    void show(){
      for (long int  i = this->deg ;i >= 0; i--)
      {this->coeff[i].show(); 
       cout << "x^" << i;
       if(i>0) {cout<<"+";}
        
      }
      cout << endl;
    }
    FiveTuple operator= (const FiveTuple &poly_b){
      this->set_to_zeros();
      for (long int  i = 0; i <= max(poly_b.deg,this->deg); i++)
      {
        this->coeff[i].numerator = poly_b.coeff[i].numerator;
        this->coeff[i].denominator = poly_b.coeff[i].denominator;
      }
      this->deg = poly_b.deg;
      return *this;

  } 
  
    
    FiveTuple operator+ (const FiveTuple &poly_b){
      FiveTuple tmp(0);  // Cn we use emoty argument in constructor?
      for (int i = 0; i <= max(poly_b.deg,this->deg); i++)
      {
        tmp.coeff[i] = this->coeff[i] + poly_b.coeff[i];
      }
      tmp.degree();
      return tmp;
    }

    FiveTuple operator- (const FiveTuple &poly_b){
      FiveTuple tmp(0);  // Cn we use emoty argument in constructor?
      for (int i = 0; i <= max(poly_b.deg,this->deg); i++)
      {
        this->coeff[i] = this->coeff[i] - poly_b.coeff[i];
      }
      
      tmp.degree();
      return tmp;
    }
    FiveTuple operator*(const rational &r){
      FiveTuple tmp(0);
      for (int i = 0; i <= this->deg; i++)
      {
        tmp.coeff[i] = this->coeff[i] *r;
        
      }
      
      return tmp;
    }
    FiveTuple operator%  (FiveTuple &poly_b){

      FiveTuple tmp(0);
      tmp = *this;
      while (tmp.deg >=poly_b.deg)
      {rational r;
      
      

       r = tmp.coeff[tmp.deg] / poly_b[poly_b.deg];
       for (long int  i = poly_b.deg; i >=0; i--)
       { 
        rational t;
        tmp.coeff[tmp.deg - i] = tmp.coeff[tmp.deg -i] - (r*poly_b[poly_b.deg -i]);
        t =this->coeff[tmp.deg -i] - (r*poly_b[poly_b.deg -i]);
        
       }
       
       tmp.degree();
       
       
        
      }
      return tmp;
      
    }
    FiveTuple operator*  ( FiveTuple  &poly_b){    
      FiveTuple tmp(0);
      for (long int i = 0; i <= this->deg + poly_b.deg; i++)
      { 
        for (long int j = 0; j <= i; j++)
        {   tmp.coeff[i] += this->coeff[j] * poly_b[i-j];
           
           
         
        }
        
        
      }
      tmp.degree();
     // cout<<tmp.deg <<endl;
      //tmp.show();
      FiveTuple cyclo_poly(0);
      rational one(1,1);
      cyclo_poly.set(one,4);
      cyclo_poly.set(one,0);
      tmp = tmp %cyclo_poly;
      //cout<<tmp.deg <<endl;
     //cyclo_poly.show();
      //tmp.show();
      
      
      return tmp;
    }
  
     
};




int main()
{   FiveTuple q(0);
    FiveTuple p(0);
    FiveTuple f(0);
    FiveTuple r(0);

    

    rational one(1,1);
    rational i(1,2);
    rational j(2,3);
    rational zero(0,1);
    rational two(2,1);
    rational minus_one(-1,1);
    rational minus_two(-2,1);



    q.set(one , 2);
    
    p.set(one , 3);
    p.set(one,1);
    p.set(minus_two, 0);
   // cout <<"n is "<< minus_two.numerator << endl;
    //p[0].show();
    
    
    
    
    //cout << p.deg<< endl;
    r = q*q;
    r =p*p;
    
    
    p.degree();
   // cout << p.deg<< endl;
    //cout << r.deg <<endl;
   // p[0].show();
  //  cout << "p is ";
 //   p.show();
    cout << "r is ";
    r.show();
    
}