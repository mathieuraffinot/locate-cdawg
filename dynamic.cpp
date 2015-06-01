#include <iostream>
#include <type_traits>

template <typename T, typename Enable = void>
class dynamic_array;


template <typename T>
class dynamic_array<T, typename std::enable_if<std::is_pod<T>::value, void>::type> {

  static_assert(std::is_pod<T>::value, "T must be POD");
  
private:
  T * tabcont;
  int size;
  int nb_elem;
  float ratio;

public:

  // basic constructor, ratio=1.5
  dynamic_array<T, typename std::enable_if<std::is_pod<T>::value, void>::type> (int size_init) {
    tabcont= (T*) malloc (sizeof(T)*size_init);
    size=size_init;
    nb_elem=0;
    ratio=1.5;
  }
  
  // changing the ratio
  dynamic_array<T, typename std::enable_if<std::is_pod<T>::value, void>::type> (int size_init, float theratio) {
    tabcont= (T*) malloc (sizeof(T)*size_init);
    size=size_init;
    nb_elem=0;
    ratio=theratio;
  }

  // realloc the main table if necessary
  void push(T newelem) {
    if (nb_elem>=size-1) {
      T * tmpcont=(T*)realloc(tabcont,(int)(size*ratio)*sizeof(T));
      if (!tmpcont)
                {
                        fprintf(stderr, "ERROR: Couldn't realloc memory for the stack of occurrences!\n");
                        exit(-1);
                }
      else { tabcont=tmpcont; size*=ratio; }
    }   
    tabcont[nb_elem++]=newelem;  
  }

  // no more elements (but memory stays)
  void makeempty() {
    size=0;
  }

  // test if the stack is empty
  int isempty() {
    return size==0;
  }
  

  // realloc the main table if necessary, but non initialized new element
  T* push_empty() {
    if (nb_elem>=size-1) {
      T * tmpcont=(T*)realloc(tabcont,(int)(size*ratio)*sizeof(T));
      if (!tmpcont)
                {
                        fprintf(stderr, "ERROR: Couldn't realloc memory for the stack of occurrences!\n");
                        exit(-1);
                }
      else { tabcont=tmpcont; size*=ratio; }
    }   
    
    return &tabcont[nb_elem++];
  }


  // pointer to the first element
  T * point_top() {
    return &tabcont[nb_elem-1];
  }

  // element on the top, does not realease
  T top() {
    return tabcont[nb_elem-1];
  }

  // element on the top, releases, but no memory liberation
  T pop() {
    return tabcont[--nb_elem];
  }

  // prints the amount of memory occupied 
  void memory() {
    std::cout << " total memory used: " << size*sizeof(T) + 2* sizeof(int) << std::endl; 
  }

  // prints the whole table (for debugging)
  void print() {
    int k;
    
    std::cout << "Taille :" << size << " ratio " << ratio <<  " - last elem " << nb_elem << " - top " << top() << std::endl;
    for (k=0;k<nb_elem;k++) std::cout << tabcont[k] << ", ";
    std::cout << " -- " << std::endl;

  }
  

  int nbelem() {
    return size;
  }

};



// the main to test
// main() {
//   int i=3;

//   dynamic_array<int> * dyn = new dynamic_array<int>(2);
  
//   dyn->push(2);
//   dyn->push(3);
//   dyn->push(5);
//   dyn->push(5);
//   dyn->push(10);
//   dyn->push(30);

//   dyn->print();

//   std::cout << " pop :" << dyn->pop() << std::endl;

//   dyn->memory();
  
// }
