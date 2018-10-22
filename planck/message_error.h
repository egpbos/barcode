/*
 *  Class for error reporting
 *
 *  Copyright (C) 2003, 2004 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_MESSAGE_ERROR_H
#define PLANCK_MESSAGE_ERROR_H

#include <iostream>
#include <string>

class Message_error
  {
  private:
    std::string msg;

  public:
    Message_error()
      : msg (std::string("Unspecified error"))
      { std::cerr<<msg<<std::endl; }

    explicit Message_error(const std::string &message)
      : msg (message) { std::cerr<<msg<<std::endl; }

    virtual const char* what() const
      { return msg.c_str(); }

    virtual ~Message_error() {}
  };

#endif
