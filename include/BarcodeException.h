/*
 * Barcode
 * Copyright E.G.P. Bos and F.S. Kitaura
 *
 * Distributed under the terms of the MIT License.
 * The full license is in the file LICENSE, distributed with this software.
 */

#pragma once
#include <stdexcept>
#include <string>
#include <utility>

// TODO: replace by simply std::runtime_error
struct BarcodeException : std::runtime_error {
    explicit BarcodeException(const std::string& s) : std::runtime_error(s) {}
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