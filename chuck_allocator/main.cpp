#include <iostream>
#include <vector>
#include "chuck_allocator.h"

int main() {
  // check allocate and construct
  chuck_allocator<int16_t> a1;
  int16_t* a = a1.allocate(1);
  a1.construct(a, 7);
  std::cout << a[0]<< '\n';
  a1.deallocate(a, 1);

// check rebind
  decltype(a1)::rebind<int64_t>::other a2;

// check allocate more than one element
  int64_t* s = a2.allocate(2);
  a2.construct(s, 100500);
  a2.construct(s + 1, -4);
  std::cout << s[0] << ' ' << s[1] << '\n';
  a1.destroy(a);

// use in standart container
  std::vector<int, chuck_allocator<int>> vec;
  vec.push_back(111111);
  vec.push_back(222121324);
  std::cout << vec[0] << ' ' << vec[1] << '\n';

// check assigmnet
  chuck_allocator<int64_t> a3;
  std::cout << a2.get_counter() << " " << a3.get_counter() << "\n";
  a3 = a2;
  std::cout << a2.get_counter() << " " << a3.get_counter() << "\n";

// check copy
  chuck_allocator<int64_t> a4(a3);
  std::cout << a4.get_counter() << " " << a3.get_counter() << "\n";

  return 0;
}
