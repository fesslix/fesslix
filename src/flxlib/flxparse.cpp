/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2025 Wolfgang Betz
 *
 * Fesslix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fesslix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fesslix.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "flxparse.h"
#include "flxdata.h"
#include "flxobjects.h"

#include <pybind11/numpy.h>  // For NumPy support


FlxFunction* parse_py_para(const std::string& para_name, py::dict config, const bool required, const tuint NumbOfPara)
{
  if (config.contains(para_name.c_str()) == false) {
    if (required) {
      std::ostringstream ssV;
      ssV << "Key '" << para_name << "' not found in Python <dict>.";
      throw FlxException_NeglectInInteractive("parse_py_para_01", ssV.str());
    } else {
      return nullptr;
    }
  }
  return parse_function(config[para_name.c_str()], "key '"+para_name+"' in Python <dict>", NumbOfPara);
}

const bool parse_py_para_as_bool(const std::string& para_name, py::dict config, const bool required, const bool def_val)
{
  if (config.contains(para_name.c_str())) {
    try {
      return py::cast<bool>(config[para_name.c_str()]);
    } catch (const py::cast_error &e) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_bool_01", "Key '"+para_name+"' in Python <dict> cannot be cast into type 'bool'.");
    }
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_bool_02", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}

const tuint parse_py_obj_as_tuint(py::object obj, std::string descr)
{
    try {
      return py::cast<tuint>(obj);
    } catch (const py::cast_error &e) {
      throw FlxException_NeglectInInteractive("parse_py_obj_as_tuint", descr + " cannot be cast into type 'unsigned int'.");
    }
}

const tuint parse_py_para_as_tuint(const std::string& para_name, py::dict config, const bool required, const tuint def_val)
{
  if (config.contains(para_name.c_str())) {
    return parse_py_obj_as_tuint(config[para_name.c_str()],"Key '"+para_name+"' in Python <dict>");
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_tuint_02", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}

const tuint parse_py_para_as_tuintNo0(const std::string& para_name, py::dict config, const bool required, const tuint def_val)
{
  const tuint val = parse_py_para_as_tuint(para_name,config,required,def_val);
  if (val<1) {
    throw FlxException_NeglectInInteractive("parse_py_para_as_tuintNo0", "Key '" + para_name + "': Value must not be negative or zero.");
  }
  return val;
}

const tulong parse_py_para_as_tulong(const std::string& para_name, py::dict config, const bool required, const tulong def_val)
{
  if (config.contains(para_name.c_str())) {
    try {
      return py::cast<tulong>(config[para_name.c_str()]);
    } catch (const py::cast_error &e) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_tulong_01", "Key '"+para_name+"' in Python <dict> cannot be cast into type 'unsigned long'.");
    }
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_tulong_02", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}

const tdouble parse_py_para_as_float(const std::string& para_name, py::dict config, const bool required, const tdouble def_val)
{
  if (config.contains(para_name.c_str())) {
    try {
      return py::cast<tdouble>(config[para_name.c_str()]);
    } catch (const py::cast_error &e) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_float_01", "Key '"+para_name+"' in Python <dict> cannot be cast into type 'float'.");
    }
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_float_02", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}

const tdouble parse_py_para_as_floatPosNo0(const std::string& para_name, py::dict config, const bool required, const tdouble def_val)
{
  const tdouble val = parse_py_para_as_float(para_name,config,required,def_val);
  if (val<=ZERO) {
    throw FlxException_NeglectInInteractive("parse_py_para_as_floatPosNo0", "Key '" + para_name + "': Value must not be negative or zero.");
  }
  return val;
}

const flxVec parse_py_para_as_flxVec(const std::string& para_name, py::dict config, const bool required, const flxVec def_val)
{
  if (config.contains(para_name.c_str())) {
    return parse_py_obj_as_flxVec(config[para_name.c_str()],"Key '"+para_name+"' in Python <dict>");
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_flxVec", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}

py::dict parse_py_obj_as_dict(py::object obj, std::string descr)
{
  if (py::isinstance<py::dict>(obj)) {
    return obj.cast<py::dict>();
  } else {
    throw FlxException_NeglectInInteractive("parse_py_obj_as_dict", descr + " cannot be cast into type '<dict>'.");
  }
}

py::list parse_py_obj_as_list(py::object obj, std::string descr)
{
  if (py::isinstance<py::list>(obj)) {
    return obj.cast<py::list>();
  } else {
    throw FlxException_NeglectInInteractive("parse_py_obj_as_list", descr + " cannot be cast into type '<list>'.");
  }
}

flxVec parse_py_obj_as_flxVec(py::object obj, std::string descr)
{
  try {
    py::array_t<tdouble> arr = py::cast<py::array_t<tdouble>>(obj);

    // Access the input data as a raw pointer
    py::buffer_info buf_info = arr.request();
    tdouble* input_ptr = static_cast<tdouble*>(buf_info.ptr);

    // Get the size of the input array
    size_t size = buf_info.size;

    return flxVec(input_ptr,size,false,false);
  } catch (const py::cast_error &e) {
    throw FlxException_NeglectInInteractive("parse_py_obj_as_flxVec", descr + " cannot be cast into type <numpy.ndarray>.");
  }
}

std::string parse_py_obj_as_string(py::object obj, std::string descr)
{
    try {
      return py::cast<std::string>(obj);
    } catch (const py::cast_error &e) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_string_01", descr + " cannot be cast into type 'string'.");
    }
}

std::string parse_py_para_as_string(const std::string& para_name, py::dict config, const bool required, const std::string def_val)
{
  if (config.contains(para_name.c_str())) {
    return parse_py_obj_as_string(config[para_name.c_str()],"Key '"+para_name+"' in Python <dict>");
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_string_02", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}

std::string parse_str_as_word(std::string strV, const bool lowercase, const bool emptyAllow, const bool numbegallow)
{
  if (lowercase) {
    std::transform(strV.begin(), strV.end(), strV.begin(), (int(*)(int)) std::tolower);
  }
  if (strV.length() == 0) {
    if (emptyAllow) {
      return "";
    } else {
      std::ostringstream ssV_2;
      ssV_2 << "Word must not be empty.";
      throw FlxException("parse_py_para_as_word_1", ssV_2.str());
    }
  }
  tuint i_start = 0;
  if (!numbegallow) {
    if ( ReadStream::getType( strV[0] ) != ReadStream::STRING ) {
      std::ostringstream ssV_2;
      ssV_2 << "Evaluated string '" << strV << "' is not of type 'word'.";
      throw FlxException("parse_py_para_as_word_2", ssV_2.str());
    }
    i_start = 1;
  }
  for ( std::string::size_type i = i_start; i < strV.length(); ++i) {
    ReadStream::InpType iT = ReadStream::getType( strV[i] );
    if ( iT != ReadStream::STRING && iT != ReadStream::NUMERAL ) {
      std::ostringstream ssV_2;
      ssV_2 << "Evaluated string '" << strV << "' is not of type 'word'.";
      throw FlxException("parse_py_para_as_word_3", ssV_2.str());
    }
  }
  return strV;
}

std::string parse_py_para_as_word(const std::string& para_name, py::dict config, const bool required, const bool lowercase, const bool emptyAllow, const bool numbegallow, const std::string def_val)
{
  std::string strV = parse_py_para_as_string(para_name,config,required,def_val);
  return parse_str_as_word(strV, lowercase, emptyAllow, numbegallow);
}

void parse_py_para_as_word_lst(std::vector<std::string>& res, const std::string& para_name, py::dict config, const bool required, const bool lowercase, const bool emptyAllow, const bool numbegallow)
{
  // ceck if 'para_name' is a key in 'config'
    if (config.contains(para_name.c_str())==false) {
      if (required) {
        throw FlxException_NeglectInInteractive("parse_py_para_as_str_lst_01", "Key '" + para_name + "' not found in Python <dict>.");
      } else {
        return;
      }
    }
  // transform to py::list
    if (py::isinstance<py::list>(config[para_name.c_str()])) {
        py::list lst = config[para_name.c_str()].cast<py::list>();
        res.reserve(res.size()+lst.size());
        for (pybind11::ssize_t i = 0; i < lst.size(); ++i) {
          py::object obj = lst[i];
          const std::string entry = parse_str_as_word(parse_py_obj_as_string(obj, "list entry"), lowercase, emptyAllow, numbegallow);
          res.push_back(entry);
        }
    } else {
      throw FlxException_NeglectInInteractive("parse_py_para_as_str_lst_02", "Key '"+para_name+"' in Python <dict> cannot be cast into type 'list'.");
    }
}

std::vector<std::string> parse_strseq_as_vec(const std::string& strseq, const char sep)
{
  std::vector<std::string> res;
  std::size_t posLast = 0;
  std::size_t pos;
  do {
    pos = strseq.find(',',posLast);
    // extract the set
      std::string sn = strseq.substr(posLast,pos-posLast);
      trim(sn);
    // add it to the list
      res.push_back(sn);
    // get the next position
      posLast = pos + 1;
  } while ( pos != std::string::npos );
  return res;
}

FlxFunction* parse_function(const std::string& funStr)
{
  return get_ReadManager()->parse_function(funStr);
}

FlxFunction* parse_function(py::object pyobj, std::string descr, const tuint NumbOfPara)
{
  return get_ReadManager()->parse_function(pyobj, descr, NumbOfPara);
}

FlxObjBase* parse_code(py::object pyobj, std::string descr)
{
  try {
    const std::string val = py::cast<std::string>(pyobj);
    return get_ReadManager()->parse_code(val);
  } catch (const py::cast_error &e) {
    throw FlxException_NeglectInInteractive("parse_code_01", descr + " cannot be cast into type 'float'.");
  }
}

template <typename T>
py::array_t<T> py_wrap_array_no_ownership(T* ptr, size_t N)
{
    // Dummy capsule — ensures py::array holds on to a base object,
    // but does not delete ptr
    py::capsule dummy_capsule(ptr, [](void* /*f*/) {
        // No deletion, memory is managed elsewhere
    });

    return py::array_t<T>(
        {N},             // shape
        {sizeof(T)}, // stride (in bytes)
        ptr,             // raw data pointer
        dummy_capsule    // base object
    );
}
template py::array_t<tfloat> py_wrap_array_no_ownership<tfloat>(tfloat*, size_t);
template py::array_t<tdouble> py_wrap_array_no_ownership<tdouble>(tdouble*, size_t);
template py::array_t<const tdouble> py_wrap_array_no_ownership<const tdouble>(const tdouble*, size_t);

py::array_t<tdouble> convert_flxVec_to_pyArray(flxVec& vec)
{
    const tuint N = vec.get_N();
    // Allocate memory for the return array
      auto res_buf = py::array_t<tdouble>(N);
    // Get the buffer info to access the underlying return data
      py::buffer_info res_buf_info = res_buf.request();
      tdouble* res_ptr = static_cast<tdouble*>(res_buf_info.ptr);
      flxVec tmp(res_ptr,N);
    // copy values
      tmp = vec;
    // Return the array
    return res_buf;
}

