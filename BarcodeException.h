/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once
#include <exception>
#include <string>

class BarcodeException : public std::exception {
  private:
    std::string s;
  public:
    BarcodeException(std::string ss) : s(ss) {}
    ~BarcodeException() throw () {}
    const char* what() const throw() { return this->s.c_str(); }
};

/* Usage example:
void Foo::Bar(){
  if(!QueryPerformanceTimer(&m_baz)){
    throw BarcodeException("it's the end of the world!");
  }
}

void Foo::Caller(){
  try{
    this->Bar();// should throw
  }catch(BarcodeException& caught){
    std::cout<<"Got "<<caught.what()<<std::endl;
  }
}
*/