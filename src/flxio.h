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

#include "flxVec.h"

#include <algorithm>
#include <stack>
#include <vector>
#include <ostream>

std::ofstream* open_stream(std::string filename);

const std::string flx_date2str(time_t date);
const std::string flx_time2str(time_t timevar);

template<class T, class U>
struct flx_is_same {
    enum { value = 0 };
};

template<class T>
struct flx_is_same<T, T> {
    enum { value = 1 };
};


//------------------------------------------------------------

class istream_warper {
  private:
    /**
    * @brief The file stream
    */
    std::istream *theStream;
    /**
    * The name of the file which is streamed.
    */
    const std::string filename;
    std::stack<int> mystack;
  public:
    istream_warper(std::istream *theStream, const std::string filename, const bool errSerious);
    ~istream_warper();
    const std::string& get_FileName() const { return filename; }
    const int get();
    const int peek();
    const bool eof();
    void putback(const int ch);
    const std::string get_line(const char delim=10);
};

//------------------------------------------------------------

struct FlxReaderPos {
  int lineNumb;
  int charNumb;
  std::string fileName;
};

//------------------------------------------------------------

const std::string ReadStream_String_End = "\n            ";
/**
* @brief A class for extracting information out of a file.
*/
class ReadStream {
  public:
    /**
    * @brief Labeles the type of a character.
    *
    * STRING: 'A'-'Z', 'a'-'z', '_'
    * NUMERAL: '0'-'9'
    * BRAKET: '{', '}', '[', ']', '(', ')'
    * QUOTE: '"', '''
    * COMMENT: '#'
    * NONE: tab, newline, blank
    * ENDOFFILE: eof
    * OTHER: any other character between '!'-'~'
    * UNKNOWN: any other character
    */
    enum InpType { STRING, NUMERAL, OTHER, BRAKET, QUOTE, NONE, COMMENT, ENDOFFILE, UNKNOWN };
    static InpType getType(char ch);
  protected:
    /**
    * @brief The file stream (warper)
    */
    istream_warper *theStream;
    /**
    * @brief How many blanks fit into a tab?
    */
    const int TabWidth;        // tab width (for calculation of column number)
    /**
    * @brief The current line number.
    */
    int lineNumb;        // current line number
    /**
    * @brief The current column number.
    */
    int charNumb;        // current column number
    /**
    * @brief if true: output will be written to pre-log
    */
    bool do_log;
    /**
    * @brief gets the next character in the stream --> and write to pre-log
    */
    void log_char(char c) { if (do_log) GlobalVar.prelog_append(c); }
    /**
    * @brief gets the next character in the stream --> and write to pre-log
    */
    const char stream_getnext_log();
  public:
    /**
    * @brief Go to the next character in the stream which is ( !( NONE || COMMENT) )
    * @returns skipped number of empty lines (without comments)
    */
    const tuint setNext ( const bool DOlog=true );
    /**
    * @brief Create a ReadStream class.
    * @param FileName The name of the file to open.
    * @param tabWidth = 8; How many blanks fit into a tab?
    * @return ReadStream
    */
    ReadStream ( const char* FileName, bool do_log=false, int tabWidth = 8, const bool errSerious=true );
    /**
    * @brief Create a ReadStream class.
    * @param strV The string to read.
    * @param tabWidth = 8; How many blanks fit into a tab?
    * @return ReadStream
    */
    ReadStream ( std::string strV, bool do_log=false, int tabWidth = 8 );
    
//     ReadStream ( std::istream& istr, int tabWidth = 8 );  // TODO: does not work correctly
    
    virtual ~ReadStream ( );
    
    const std::string& get_FileName() { return theStream->get_FileName(); }
    /**
    * @brief Indicates the position where an error occurred.
    * @return String describing the current position in text file (line number, column number, name of the file).
    */
    const std::string getCurrentPos ();
    void getCurrentPos(FlxReaderPos& pos);
    static const std::string write_ReaderPos(const FlxReaderPos& pos);
    /**
    * @brief Returns the type of the next character in the stream.
    * @return InpType
    */
    const InpType getNextType ();

    /**
    * @brief Returns the next character in the string.
    * @param DOsetNext = false; If this parameter is true, 'setNext()' is executed afterwards.
    * @return char
    * @see setNext()
    */
    const char getChar ( const bool DOsetNext = true, const bool DOlog=true );
    
    /**
    * @brief Returns the next character in the string - and throws an error if it is not the expected one.
    * @param theChar The expected character.
    * @return char
    * @throw FlxException
    */
    const char getChar ( const char theChar, const bool errSerious = true, const bool DOsetNext = true );

    /**
    * @brief Returns the next character in the stream without moving the current position in the stream.
    * @return char
    */
    const char whatIsNextChar ( ) { return theStream->peek(); }
    /**
    * @brief Returns the next string of length 'length' in the stream without moving the current position in the stream.
    * allowed are a-z A-z 0-1 _
    * @return char
    */
    const std::string whatIsNextString ( const int length, const bool lowercase = false );
    
    /**
    * @brief returns the next line (it acutally reads it)!!! ,,, name is missleading
    */
    const std::string getNextLine( const bool DOsetNext = true );

    /**
    * @brief Returns the next word in the string.
    * 
    * The fist character of the word has to be of the type STRING.
    * For all other character the types STRING and NUMERAL are valid.
    * A character of any other type indicates the end of the word.
    * The function 'setNext()' is executed afterwards.
    * @return A string containing the extracted word.
    * @see setNext()
    * @throw FlxException
    */
    const std::string getWord (const bool lowercase, const bool errSerious = true, const bool numbeg=false, const bool dblpnt=false);
    /**
    * @brief Checks if the next word in the stream is 'strWord'
    */
    const std::string getWord(const char* strWordC, const bool errSerious = true);
    /**
    * @brief Checks if the next expression in the stream is 'strExpr'
    */
    void getExpr(const char* strWordC, const bool errSerious = true);

    /**
    * @brief Returns the next bool-value in the string.
    * 
    * Allowed is:
    *  'true' or 'false'
    *  '1' or '0'
    * @return true or false
    * @see setNext()
    * @throw FlxException
    */
    const bool getBool(const bool errSerious = true);

    /**
    * @brief Returns the next word in the string.
    * 
    * The text TEXT must have the following format: "TEXT"
    * @return A string containing the extracted text.
    * @see setNext()
    * @throw FlxException
    */
    const std::string getText (const bool MultiLine = false, const bool errSerious = true);
    
    /**
    * @brief Checks if the next character can be the starting character of a number.
    * @return true or false
    */
    const bool nextCanBeNumber();
    
    /**
    * @brief forgets everything before strexpr & the expression itself
    * @return true if found, false if an error occurs
    */
    const bool set_after_expr(const std::string& strexpr, const bool DOsetNext = true);
    /** 
    * @brief ignore all characters until the first one matches type
    * @returns false if EOF
    */
    const bool ignore_until(const InpType typ);
    /** 
    * @brief ignore all characters until the closing bracket on same level
    * @param bracket ')', '}', ']'
    */
    void ignore_bracket(const char bracket);

    time_t getDate(const bool errSerious = true);
    time_t getTime(struct tm date_tm, const bool errSerious = true);
    
    /**
    * @brief Returns the next number in the string.
    * 
    * The format of the string describing the number must be of this format:
    * ['+','-']'0'-'9'...'0'-'9'[.'0'-'9'...'0'-'9'][('e', 'E')['+','-']'0'-'9'...'0'-'9']
    * Note: the use of 'e' or 'E' to extract an integer value can cause problems.
    * The function 'setNext()' is executed afterwards.
    * @return The extracted numeric value.
    * @see setNext()
    * @throw FlxException
    */
    const tdouble get_Double(const bool errSerious = true);
    /** 
    * @returns a floating-point number that is given in the format 9.487,61
    */
    const tdouble get_GermanFloat(const bool errSerious = true);
    
    /** 
    * @returns an unsigned integer number
    */
    template <class T> const T get_UInt( const bool errSerious=true, const bool doSetNext=true ) {
        bool read_number = false;
        // read the number before the dot
          T T1(0);        // the current number retrieved
          char ch = theStream->peek();
          while (ch>=48&&ch<=57) {
            T1 *= 10;
            T1 += T(ch-48);
            read_number = true;
            getChar(false);
            ch = theStream->peek();
          }
        // check for errors
          if (!read_number) {
            std::ostringstream ssV;
            ssV << "A number is required at this point.";
            FlxError(errSerious,"ReadStream::getNumber", ssV.str(), getCurrentPos());
          }
        // set reader to next ...
          if (doSetNext) setNext();
        return T1;
    }

};

typedef ReadStream* ReadStream_ptr;



//------------------------------------------------------------


class FlxReaderBase {
  protected:
    static ReadStream *reader;
  public:
    virtual ~FlxReaderBase() {};
};


//------------------------------------------------------------



class FlxIstream {
  protected:
    /**
    * @brief the name of the inputReader
    */
    std::string name;
    /**
    * @brief generate an error if end of file is not respected (Warning by default)
    */
    bool errEOF;
    
    
    /**
    * @brief this function is executed if the end of the file is reached
    */
    void reachedEOF();
  public:
    FlxIstream(std::string nameV, bool errEOFv);
    virtual ~FlxIstream() {}
    
    virtual void copyStream(FlxIstream* rhsB, bool errSerious=true, const bool err_check=true) = 0;
    /**
    * @brief returns the next value
    * @return false if EOF is reached
    */
    virtual bool get_value(tdouble& v, bool suppressErr = false) = 0;
    /**
    * @brief returns the next value
    * @param v the vector to assemble
    * @param lastIndex if EOF occurred, it contains the index where the error occurred within the vector
    * @return false if EOF is reached
    */
    const bool get_vec(flxVec& v, tuint &lastIndex, const bool suppressErr = false);
};

class FlxIstream_file : public FlxIstream {
  protected:
    /**
    * @brief random sequence input: reader
    */
    ReadStream* iReader;
    /**
    * @brief random sequence input: blocksize
    */
    tuint blocksize;
    /**
    * @brief random sequence input: sequence vector
    */
    tVec* SeqVec;
    /**
    * @brief current index in SeqVec
    */
    tuint index;
    /**
    * @brief last valid index in SeqVec
    */
    tuint lastindex;
    /**
    * @brief sets the input to the next number (ignoring ',' and ';')
    */
    void set_next();
    /**
    * @brief read another block
    */
    virtual void read_block();
    /**
    * @brief check for end of file
    */
    bool check_eof(ReadStream_ptr& rp);
    
    // reading specific columns
      /**
      * @brief the number of columns in the file
      */
      tuint Cnumb;
      /**
      * @brief the column numbers to read - starts counting with 1 !!!
      */
      std::vector<tuint> Cvec;
      /**
      * @brief the current column-number of the input (the next to read)
      */
      tuint curCol;
      /**
      * @brief the current index in the vector Cvec; (the next to check)
      */
      tuint curVecIdx;
      
  public:
    FlxIstream_file(std::string nameV, ReadStream* iReaderV, bool errEOFv, tuint blocksizeV, const tuint Cnumb, const std::vector<tuint> &Cvec, const bool do_read_block=true);
    virtual ~FlxIstream_file();
   
    virtual void copyStream(FlxIstream* rhsB, bool errSerious=true, const bool err_check=true);
    virtual bool get_value(tdouble& v, bool suppressErr = false);
};

class FlxIstream_file_binary : public FlxIstream_file {
  private:
    std::ifstream* file;
    std::streampos size;
    bool binaryfloat;                // float or double values?
    /**
    * @brief read another block
    */
    virtual void read_block();
    
  public:
    FlxIstream_file_binary(std::string nameV, const std::string& fileName, bool errEOFv, tuint blocksizeV, const tuint Cnumb, const std::vector<tuint> &Cvec, const bool binaryfloat);
    ~FlxIstream_file_binary();
   
    virtual void copyStream(FlxIstream* rhsB, bool errSerious=true, const bool err_check=true);
    
    const size_t get_N_numbers();
};

class FlxIstream_file_combine : public FlxIstream_file {
  private:
    std::vector<ReadStream*> iReader_vec;
    flxVec* weight_vecP;
    /**
    * @brief read another block
    */
    virtual void read_block();
    
  public:
    FlxIstream_file_combine(std::string nameV, std::vector<ReadStream*>& iReader_vec, flxVec& weight_vec, bool errEOFv, tuint blocksizeV, const tuint Cnumb, const std::vector<tuint> &Cvec);
    ~FlxIstream_file_combine();
   
    virtual void copyStream(FlxIstream* rhsB, bool errSerious=true, const bool err_check=true);
};

class FlxIstream_vector : public FlxIstream {
  private:
    std::vector<tdouble> iVec;
    /**
    * @brief current index in iVec
    */
    tulong index;
    /**
    * @brief number of entries in iVec
    */
    tulong numbEl;
    
  public:
    FlxIstream_vector(std::string nameV, FlxIstream* is, bool errEOFv, const tulong Nreserve);
    
    void copyStream(FlxIstream* rhsB, bool errSerious=true, const bool err_check=true);
    bool get_value(tdouble& v, bool suppressErr = false);
    
    /**
    * @brief append a number to the end of the stream
    */
    void appendNumber(const tdouble &d);
    /**
    * @brief delete all numbers in the stream
    */
    void clear();
    
    /**
    * @brief sorts the vector and resets the stream
    */
    void sortStream();
    /**
    * @brief resets the stream
    */
    void reset_stream() { index = 0; }
    /**
    * @brief returns the total number of entries in the stream
    */
    tulong get_total_size() { return numbEl; }
    /**
    * @brief returns the first value in the set
    */
    tdouble get_first_value() { return iVec.empty()?ZERO:iVec[0]; }
    /**
    * @brief returns the last value in the set
    */
    tdouble get_last_value() { return iVec.empty()?ZERO:iVec[iVec.size()-1]; }
    
    const tdouble* get_tmpPtr() { return &(iVec[0]); }
};


//------------------------------------------------------------

class flxBufAlloc : public std::streambuf
{
  private:
    typedef std::streambuf base_type;
    typedef base_type::traits_type traits_type;

    flxBufAlloc( const flxBufAlloc& );
    flxBufAlloc& operator=( const flxBufAlloc& ); // Kopieren verhindern

    const ostreamp& os1;
    const ostreamp& os2;

  protected:
    virtual int_type overflow( int_type c = traits_type::eof() )  {   
        std::streambuf* m_sb1 = os1->rdbuf();
        std::streambuf* m_sb2 = os2->rdbuf();
        if( m_sb1 && m_sb2 
            && !traits_type::eq_int_type( m_sb1->sputc( c ), traits_type::eof() )
            && !traits_type::eq_int_type( m_sb2->sputc( c ), traits_type::eof() )
            )
            return traits_type::not_eof( c );
        return traits_type::eof();  // einer der Streambufs tut nicht
    }

  public:
    flxBufAlloc ( const ostreamp& os1, const ostreamp& os2 ) : os1(os1), os2(os2) { }

};

class flxDummyAlloc : public std::streambuf
{
  private:
    typedef std::streambuf base_type;
    typedef base_type::traits_type traits_type;

    flxDummyAlloc( const flxDummyAlloc& );
    flxDummyAlloc& operator=( const flxDummyAlloc& ); // Kopieren verhindern

  protected:
    virtual int_type overflow( int_type c = traits_type::eof() )  {   
       return c;
    }

  public:
    flxDummyAlloc ( ) { }

};

class flxStreamAlloc : public std::ostream {
  private:
    flxBufAlloc theBuf;

  public:
    flxStreamAlloc ( const ostreamp& os1, const ostreamp& os2 ) : std::ostream(&theBuf), theBuf(os1, os2) { }
};

class flxDummyOstream : public std::ostream {
  private:
    flxDummyAlloc theBuf;

  public:
    flxDummyOstream ( ) :std::ostream(&theBuf) {  }
};



