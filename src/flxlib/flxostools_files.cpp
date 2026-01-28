/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2026 Wolfgang Betz
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

#define FLXLIB_CPP

#include "flxostools_files.h"
#include "flxglobal.h"
#include <stdio.h>


const bool fileIsReadable(const std::string& fileName)
{
  FILE* ofile = fopen(fileName.c_str(),"r");
  const bool res = (ofile!=NULL);
  if (res) fclose(ofile);
  return res;
}


// -----------------------------------------------------------------------------------------------------------------
#if FLX_BOOST_FS
// -----------------------------------------------------------------------------------------------------------------

#include <boost/filesystem.hpp>
  namespace fs = boost::filesystem;


// taken from http://stackoverflow.com/questions/8593608/how-can-i-copy-a-directory-using-boost-filesystem_error
//   and modified
const bool copyDir( boost::filesystem::path const & source, boost::filesystem::path const & destination )
{
    try
    {
        // Check whether the function call is valid
        if( !fs::exists(source) || !fs::is_directory(source) )
        {
            GlobalVar.alert.alert("copyDir_1","Source directory " + source.string() + " does not exist or is not a directory.");
            return false;
        }
        if(fs::exists(destination))
        {
            GlobalVar.alert.alert("copyDir_2","Destination directory " + destination.string() + " already exists.");
            return false;
        }
        // Create the destination directory
        if(!fs::create_directory(destination))
        {
            GlobalVar.alert.alert("copyDir_3","Unable to create destination directory " + destination.string());
            return false;
        }
    }
    catch(fs::filesystem_error const & e)
    {
        GlobalVar.alert.alert("copyDir_4",e.what());
        return false;
    }
    // Iterate through the source directory
    for( fs::directory_iterator file(source); file != fs::directory_iterator(); ++file )
    {
        try
        {
            fs::path current(file->path());
            if(fs::is_directory(current))
            {
                // Found directory: Recursion
                if(
                    !copyDir(
                        current,
                        destination / current.filename()
                    )
                )
                {
                    return false;
                }
            }
            else
            {
                // Found file: Copy
                fs::copy_file(
                    current,
                    destination / current.filename()
                );
            }
        }
        catch(fs::filesystem_error const & e)
        {
            GlobalVar.alert.alert("copyDir_5",e.what());
        }
    }
    return true;
}


void copyDir( const std::string& source_dir, const std::string& dest_dir )
{
  fs::path const & source = source_dir;
  fs::path const & destination = dest_dir;
  
  if (copyDir(source,destination) == false) {
    throw FlxException("copyDir_100","The directory '"+source_dir+"' could not be copied to '"+dest_dir+"'.");
  }
}

void copyFile(const std::string& source, const std::string& dest, const bool overwrite)
{
  try {
    if (overwrite) {
      fs::copy_file(source,dest,boost::filesystem::copy_options::overwrite_existing);
    } else {
      fs::copy_file(source,dest);
    }
  } catch(fs::filesystem_error const & e) {
    throw FlxException("copyFile_100","File '"+source+"' could not be copied to '"+dest+"'.",e.what());
  }
}

void moveFile(const std::string& source, const std::string& dest)
{
  try {
    fs::rename(source,dest);
  } catch(fs::filesystem_error const & e) {
    throw FlxException("moveFile_100","File '"+source+"' could not be moved to '"+dest+"'.",e.what());
  }
}

const tuint delDir(const std::string& dir)
{
  try {
    return fs::remove_all( dir );
  } catch(fs::filesystem_error const & e) {
    throw FlxException("delDir_100","Directory '"+dir+"' could not be removed.",e.what());
  }
}

const bool existsDir(const std::string& dir)
{
  try {
    return fs::exists(dir) && fs::is_directory(dir);
  } catch(fs::filesystem_error const & e) {
    throw FlxException("existsDir_100",e.what());
  }
}

const bool createDir(const std::string& dir)
{
  try {
    return fs::create_directory(dir);
  } catch(fs::filesystem_error const & e) {
    throw FlxException("createDir_100",e.what());
  }
}

// -----------------------------------------------------------------------------------------------------------------
#if FLX_BOOST_RX
// -----------------------------------------------------------------------------------------------------------------

#include <boost/regex.hpp>

void getFiles( const std::string& sourcePath, const std::string& sourcePattern, std::vector< std::string >& match_files ) {
  try {
    fs::directory_iterator end_itr; 
    for( fs::directory_iterator i( sourcePath ); i != end_itr; ++i )
    {
        // Skip if not a file
        if( !fs::is_regular_file( i->status() ) ) continue;

        boost::smatch what;
        const boost::regex my_filter( sourcePattern );
        const std::string cfn = i->path().filename().string();
        // Skip if no match
        if( !boost::regex_match(cfn, what, my_filter ) ) continue;

        // File matches, store it
        match_files.push_back( cfn );
    }
  } catch (boost::bad_expression &e) {
    throw FlxException("getFiles",e.what());
  }
}


// -----------------------------------------------------------------------------------------------------------------
#else   // FLX_BOOST_RX
// -----------------------------------------------------------------------------------------------------------------
  
void getFiles( const std::string& sourcePath, const std::string& sourcePattern, std::vector< std::string >& match_files ) {
  throw FlxException("getFiles_300","Method not available.");
}

// -----------------------------------------------------------------------------------------------------------------
#endif  // FLX_BOOST_RX
// -----------------------------------------------------------------------------------------------------------------


// -----------------------------------------------------------------------------------------------------------------
#else   // FLX_BOOST_FS
// -----------------------------------------------------------------------------------------------------------------

void copyDir( const std::string& source_dir, const std::string& dest_dir, const bool overwrite )
{
  throw FlxException("copyDir_200","Method not available.");
}

void copyFile(const std::string& source, const std::string& dest, const bool overwrite)
{
  throw FlxException("copyFile_200","Method not available.");
}

void moveFile(const std::string& source, const std::string& dest)
{
  throw FlxException("moveFile_200","Method not available.");
}

const tuint delDir(const std::string& dir)
{
  throw FlxException("delDir_200","Method not available.");
}

const bool existsDir(const std::string& dir)
{
  throw FlxException("existsDir_200","Method not available.");
}

const bool createDir(const std::string& dir)
{
  throw FlxException("createDir_200","Method not available.");
}

void getFiles( const std::string& sourcePath, const std::string& sourcePattern, std::vector< std::string >& match_files ) {
  throw FlxException("getFiles_200","Method not available.");
}

// -----------------------------------------------------------------------------------------------------------------
#endif  // FLX_BOOST_FS
// -----------------------------------------------------------------------------------------------------------------
