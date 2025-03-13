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

#pragma once

#include <string>
#include <sstream>
#include <ostream>

#if defined(_MSC_VER)
  #include "flxlib_export.h"
#else
  #define FLXLIB_EXPORT
#endif

typedef std::ostream* ostreamp;


/**
* @brief A class for exception handling.
*/
class FLXLIB_EXPORT FlxException : public std::exception {
  private:
    std::string errNumber;
    std::string titel;
    std::string msg;
    std::string full_msg;
  public:
    /**
    * @brief Create a FlxException class.
    * @param errnumber Where did the error occur? (Format: Class::Function_InternNumber)
    * @param Titel A short message describing the error.
    * @param Msg More error information.
    * @return FlxException
    */
    FlxException( std::string errnumber, std::string Titel = "", std::string Msg = "" );
    /**
    * @brief Get error information (usefull for catch-blocks).
    * @return string describing the error.
    */
    const char* what() const noexcept override;
};

/**
* @brief A class for exception handling.
*/
class FLXLIB_EXPORT FlxException_NeglectInInteractive : public FlxException {
  public:
    /**
    * @brief Create a FlxException_NeglectInInteractive class.
    * @param errnumber Where did the error occur? (Format: Class::Function_InternNumber)
    * @param Titel A short message describing the error.
    * @param Msg More error information.
    * @return FlxException
    */
    FlxException_NeglectInInteractive ( std::string errnumber, std::string Titel = "", std::string Msg = "" ) 
      : FlxException(errnumber, Titel, Msg) { }
};

/**
* @brief A class for exception handling.
*/
class FLXLIB_EXPORT FlxException_NotImplemented : public FlxException {
  public:
    /**
    * @brief Create a FlxException_NotImplemented class.
    * @param errnumber Where did the error occur? (Format: Class::Function_InternNumber).
    * @return FlxException
    */
    FlxException_NotImplemented ( std::string errnumber ) 
      : FlxException(errnumber, "Feature not implemented", "The requested feature has not yet been implemented.") { }
};

/**
* @brief A class for exception handling.
*/
class FLXLIB_EXPORT FlxException_Crude : public FlxException {
  public:
    /**
    * @brief Create a FlxException_Crude class.
    * @param errnumber Where did the error occur? (Format: Class::Function_InternNumber).
    * @return FlxException
    */
    FlxException_Crude ( std::string errnumber ) 
      : FlxException(errnumber, "ERROR", "Actually, this error should have never occurred ...") { }
};

/**
* @brief A class for exception handling (internal errors - to be catched)
*/
class FLXLIB_EXPORT FlxException_math : public FlxException {
  public:
    /**
    * @brief Create a FlxException_math class.
    * @param errnumber Where did the error occur? (Format: Class::Function_InternNumber).
    * @return FlxException
    */
    FlxException_math ( std::string errnumber, std::string Titel = "", std::string Msg = "" )
      : FlxException(errnumber, Titel, Msg) { }
};

/**
* @brief Used to terminate the program
*/
class FlxEndE {
  public:
    FlxEndE ( ) { }
};

/**
* @brief Used to exit from a certain read level
*/
class FlxExitE {
  public:
    FlxExitE ( ) { }
};

class FlxReturnBreakContinue_baseE : public FlxException {  // intended as common basis for 'break', 'continue' and 'return' statements
  public:
    FlxReturnBreakContinue_baseE ( std::string typestr ) : FlxException(typestr, "'"+typestr+"' executed by user outside of designated environment.") { }
};

/**
* @brief Used to return from a certain code-block
*/
class FlxReturnE : public FlxReturnBreakContinue_baseE {
  public:
    FlxReturnE ( ) : FlxReturnBreakContinue_baseE("return") { }
};

/**
* @brief Used to return from a certain code-block
*/
class FlxContinueE : public FlxReturnBreakContinue_baseE {
  public:
    FlxContinueE ( ) : FlxReturnBreakContinue_baseE("continue") { }
};

/**
* @brief Used to return from a certain code-block
*/
class FlxBreakE : public FlxReturnBreakContinue_baseE {
  public:
    FlxBreakE ( ) : FlxReturnBreakContinue_baseE("break") { }
};


/**
* @brief A class for printing alerts.
*/
class FlxAlert {
  private:
    static const ostreamp* cerr;
  public:
    /**
    * @brief Create a FlxAlert class.
    * @param &cerrV output_error_stream (Format: ostream).
    */
    FlxAlert () {}
    static void set_err_stream(const ostreamp& cerrV) { cerr = &cerrV; }
//     static const ostreamp& get_err_stream() { return *cerr; }
    /**
    * @brief Output an alert.
    * @param alertID Where did the error occur? (Format: Class::Function_InternNumber)
    * @param alertStr A short message describing the alert.
    */
    void alert(const std::string& alertID, const std::string& alertStr) const;
};


