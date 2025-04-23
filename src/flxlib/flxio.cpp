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


#include "flxio.h"

#include <fstream>
#include <cstring>
#include <cctype>
#include <boost/concept_check.hpp>


using namespace std;

ReadStream* FlxReaderBase::reader = NULL;


ofstream* open_stream(string filename)
{
  std::ofstream *theStream = new std::ofstream(filename.c_str());
  if ( ! theStream->is_open() || theStream == NULL ) {
    std::ostringstream ssV;
    ssV << "File (" << filename << ") could not be opened.";
    throw FlxException("flxio.h::open_stream", ssV.str() );
  }
  return theStream;
}

const std::string flx_date2str(time_t date)
{
  std::ostringstream ssV;
  struct tm * cdate;
  cdate = localtime(&(date));
  ssV << cdate->tm_year+1900 << '-' << cdate->tm_mon+1 << '-' << cdate->tm_mday;
  return ssV.str();
}

const string flx_time2str(time_t timevar)
{
  std::ostringstream ssV;
  struct tm * ctime;
  ctime = localtime(&(timevar));
  ssV << ctime->tm_hour << ':' << ctime->tm_min << ':' << ctime->tm_sec;
  return ssV.str();
}

istream_warper::istream_warper(istream* theStream, const string filename, const bool errSerious): theStream(theStream),filename(filename)
{
  if (filename=="(String)") return;
  std::ifstream *filstr = dynamic_cast<std::ifstream*> (theStream);
  if ( ! filstr->is_open() || filstr == NULL ) {
    ostringstream ssV;
    ssV << "File (" << filename << ") could not be opened.";
    // free memory
      std::string str1 = filename;
      if (str1 != "(cin)") delete theStream; 
    FlxError(errSerious,"istream_warper::istream_warper_1", ssV.str());
  }
}

istream_warper::~istream_warper()
{
  std::ifstream *filstr = dynamic_cast<std::ifstream*> (theStream);
  if ( filstr != NULL) {
    filstr->close();
  }
  const std::string str1 = filename;
  if (str1 != "(cin)") {
    if (theStream!=NULL) delete theStream; 
  }
}

const int istream_warper::get()
{
  if (!mystack.empty()) {
    char ch = mystack.top();
    mystack.pop();
    return ch;
  } else {
    return theStream->get();
  }
}

const string istream_warper::get_line(const char delim)
{
  std::string st;
  char ch;
  while ((ch=get())!=delim) {
    if (theStream->eof()) return st;
    st += ch;
  }
  return st;
}

const int istream_warper::peek()
{
  if (!mystack.empty()) {
    return mystack.top();
  } else {
    return theStream->peek();
  }
}

const bool istream_warper::eof()
{
  if (!mystack.empty()) {
    return false;
  } else {
    return theStream->eof();
  }
}

void istream_warper::putback(const int ch)
{
  if (ch < 0) {
    ostringstream ssV;
    ssV << "ERROR (" << ch << ")" ;
    throw FlxException("istream_warper::putback", ssV.str());
  }
  mystack.push(ch);
}

ReadStream::ReadStream ( const char* FileName, bool do_log, int tabWidth, const bool errSerious ) 
: theStream(new istream_warper(new std::ifstream(FileName, std::ios_base::in | std::ios_base::binary ), FileName,errSerious)), TabWidth(tabWidth), lineNumb(1), charNumb(0),do_log(do_log) { 
  setNext(); 
}

ReadStream::ReadStream ( string strV, bool do_log, int tabWidth ) 
: theStream(new istream_warper(new std::istringstream(strV+ReadStream_String_End),"(String)",true)), TabWidth(tabWidth), lineNumb(1), charNumb(0),do_log(do_log) { 
  setNext(); 
}

// ReadStream::ReadStream(istream& istr, int tabWidth)
// : theStream(&istr), filename("(cin)"), TabWidth(tabWidth),lineNumb(1),charNumb(0)
// {
//   setNext();
// }


ReadStream::~ReadStream ( ) { 
  delete theStream;
}

ReadStream::InpType ReadStream::getType(char ch)
{
  if ( int(ch) < 0 ) return ENDOFFILE;
  switch (ch) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
      break;
    case 9:                // tab
    case 10:                // newline
      return NONE;
    case 11:
    case 12:
      break;
    case 13:                // carriage return
      return NONE;
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
    case 19:
    case 20:
    case 21:
    case 22:
    case 23:
    case 24:
    case 25:
    case 26:
    case 27:
    case 28:
    case 29:
    case 30:
    case 31:
      break;
    case 32:                // blank
      return NONE;
    case 33:
      return OTHER;
    case 34:
      return QUOTE;
    case 35:
      return COMMENT;
    case 36:
    case 37:
    case 38:
    case 39:
      return OTHER;
    case 40:
    case 41:
      return BRAKET;
    case 42:
    case 43:
    case 44:
    case 45:
    case 46:
    case 47:
      return OTHER;
    case 48:
    case 49:
    case 50:
    case 51:
    case 52:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57:
      return NUMERAL;
    case 58:
    case 59:
    case 60:
    case 61:
    case 62:
    case 63:
    case 64:
      return OTHER;
    case 65:
    case 66:
    case 67:
    case 68:
    case 69:
    case 70:
    case 71:
    case 72:
    case 73:
    case 74:
    case 75:
    case 76:
    case 77:
    case 78:
    case 79:
    case 80:
    case 81:
    case 82:
    case 83:
    case 84:
    case 85:
    case 86:
    case 87:
    case 88:
    case 89:
    case 90:
      return STRING;
    case 91:
      return BRAKET;
    case 92:
      return OTHER;
    case 93:
      return BRAKET;
    case 94:
      return OTHER;
    case 95:
      return STRING;
    case 96:
      return QUOTE;
    case 97:
    case 98:
    case 99:
    case 100:
    case 101:
    case 102:
    case 103:
    case 104:
    case 105:
    case 106:
    case 107:
    case 108:
    case 109:
    case 110:
    case 111:
    case 112:
    case 113:
    case 114:
    case 115:
    case 116:
    case 117:
    case 118:
    case 119:
    case 120:
    case 121:
    case 122:
      return STRING;
    case 123:
      return BRAKET;
    case 124:
      return OTHER;
    case 125:
      return BRAKET;
    case 126:
      return OTHER;
    case 127:
      return UNKNOWN;
    default:
      return UNKNOWN;
  }
  return UNKNOWN;
}

const ReadStream::InpType ReadStream::getNextType () {
  if ( theStream->eof() ) return ENDOFFILE;
  return getType( theStream->peek() );
}

const bool ReadStream::check_eof()
{
  if (getNextType() == ReadStream::ENDOFFILE) {
    return true;
  }
  else return false;
}

const std::string ReadStream::getCurrentPos () {
  ostringstream ssV;
  ssV << "Line: " << lineNumb << "; Column: " << charNumb << "; File: " << theStream->get_FileName();
  return ssV.str();
}

void ReadStream::getCurrentPos(FlxReaderPos& pos)
{
  pos.lineNumb = lineNumb;
  pos.charNumb = charNumb;
  pos.fileName = theStream->get_FileName();
}

const std::string ReadStream::write_ReaderPos(const FlxReaderPos& pos)
{
  ostringstream ssV;
  ssV << "Line: " << pos.lineNumb << "; Column: " << pos.charNumb << "; File: " << pos.fileName;
  return ssV.str();
}

const char ReadStream::stream_getnext_log()
{
  const char t = theStream->get();
  log_char(t);
  return t;
}

const tuint ReadStream::setNext ( const bool DOlog ) {
  char ch;
  tuint elc = 0;        // empty line count
  InpType it = getNextType();
  while ( it == NONE || it == COMMENT ) {
    if ( it == COMMENT ) {                        // ... if there is a comment
      if (whatIsNextString(2)=="#!") {
        return elc;                // TODO: this might cause a bug!!!
      } else {
        theStream->get_line();
      }
      lineNumb++; 
      charNumb=0;
    } else {                                                 // ... if it is NOT a comment
      if ( (ch=theStream->get()) == 10 ) {                  // new line
        lineNumb++; ++elc;
        charNumb=0;
      } else if ( ch == 9 ) {                                // tab
        charNumb += TabWidth;
      } else if ( ch == 13 ) {
      } else {
        charNumb++;
      }
      if (DOlog) log_char(ch);
    }
    it = getNextType();
  }
  return elc;
}

const char ReadStream::getChar ( const bool DOsetNext, const bool DOlog ) {
  const char ch = (DOlog?stream_getnext_log():theStream->get());
  if (getType(ch)==NONE) {
    if ( ch == 10 ) {                         // newline
      lineNumb++;
      charNumb=0;
    } else if ( ch == 9) {                // tab
      charNumb += TabWidth;
    } else if ( ch == 13 ) {
    } else {
      charNumb++;
    }
  } else {
    charNumb++;
  }
  if ( DOsetNext ) setNext(DOlog);
  return ch;
}

const char ReadStream::getChar ( const char theChar, const bool errSerious, const bool DOsetNext ) {
  const char ch = getChar(DOsetNext);
  if ( ch != theChar ) {
    ostringstream ssV;
    ssV << "Expected '" << theChar << "' (and NOT '" << ch << "' [" << int(ch) << "])";
    FlxError(errSerious,"ReadStream::getChar_1", ssV.str(), getCurrentPos());
  }
  return ch;
}

const std::string ReadStream::whatIsNextString( const int length, const bool lowercase) {
  if ( length <= 0 ) {
    ostringstream ssV;
    ssV << "'length' has to be greater than zero.";
    throw FlxException("ReadStream::whatIsNextString_1", ssV.str() );
  }
  // get next string
    #ifndef FLX_CV_1
      int i1[length];
      char ch1[length+1];
    #else
      int *i1 = new int[length];
      char *ch1 = new char[length+1];
    #endif
    if (theStream->eof()) {
      return "";
    }
    for (int i = 0; i < length; ++i) {
      i1[i] = theStream->get();
      if (theStream->eof()) {
        ch1[i] = '\0';
        break;
      }
      ch1[i] = i1[i];
    }
    ch1[length] = '\0';
    string str1(ch1);
  // write it back
    for (int i = std::strlen(ch1)-1;i>=0;--i) {
      if (ch1[i] != 0) theStream->putback(i1[i]);
    }
  #ifdef FLX_CV_1
    delete [] i1;
    delete [] ch1;
  #endif
  if (lowercase) {
    std::transform(str1.begin(), str1.end(), str1.begin(), (int(*)(int)) std::tolower);
  }
  return str1;
}

const string ReadStream::getNextLine(const bool DOsetNext)
{
  lineNumb++; 
  charNumb=0;
  const std::string ts = theStream->get_line();
  if (DOsetNext) setNext();
  return ts;
}

const bool ReadStream::set_after_expr(const string& strexpr, const bool DOsetNext)
{
  const size_t length = strexpr.length();
  size_t pos = 0;
  while(pos<length) {
    // extracts the next character
      if (theStream->eof()) return false;
      char c = theStream->get();
      switch (c) {
        case 9:
          charNumb+=TabWidth;
          break;
        case 10:
          lineNumb++;
          charNumb=0;
          break;
        case 13:
          break;
        default:
          charNumb++;
      }
    if (c==strexpr[pos]) ++pos;
    else pos = 0;
  }
  if (DOsetNext) {
    setNext();
  }
  return true;
}

const bool ReadStream::ignore_until(const ReadStream::InpType typ)
{
  while(getNextType()!=typ) {
    if (theStream->eof()) return false;
    theStream->get();
  }
  return true;
}

void ReadStream::ignore_bracket(const char bracket)
{
  tuint level = 0;
  char obrack;
  if (bracket==')') {
    obrack = '(';
  } else if (bracket==']') {
    obrack = '[';
  } else if (bracket=='}') {
    obrack = '{';
  } else {
    throw FlxException_Crude("ReadStream::ignore_bracket_01");
  }
  while (true) {
    if (getNextType()==ReadStream::ENDOFFILE) {
      throw FlxException("ReadStream::ignore_bracket_02","Reach end of file.");
    }
    const char c = whatIsNextChar();
    if (c==obrack) {
      ++level;
      getChar();
    } else if (c==bracket) {
      if (level<=0) {
          getChar();
          return;
        } else {
          --level;
          getChar();
        }
    } else if (c=='"') {
      getText(true);
    } else {
      getChar();
    }
  };
}

void ReadStream::getExpr(const char* strWordC, const bool errSerious)
{
  #ifndef FLX_CV_1
    char chS[std::strlen(strWordC)];
  #else
    char *chS = new char[std::strlen(strWordC)];
  #endif
  int i = 0;
  while (strWordC[i]!='\0') {
    chS[i] = theStream->get();
    switch (chS[i]) {
      case 9:
        charNumb+=TabWidth;
        break;
      case 10:
        lineNumb++;
        charNumb=0;
        break;
      case 13:
        break;
      default:
        charNumb++;
    }
    ++i;
  }
  chS[i]='\0';
  
  std::string str1(chS);
  std::string str2(strWordC);
  
  #ifdef FLX_CV_1
    delete [] chS;
  #endif
  if (str1!=str2) {
    ostringstream ssV_2;
    ssV_2 << "Expected '" << str2 << "' and not '" << str1 << "'.";
    FlxError(errSerious,"ReadStream::getExpr", ssV_2.str(), getCurrentPos());
  }
  
  if (do_log) GlobalVar.prelog_append(str1);
}

const std::string ReadStream::getWord(const char* strWordC, const bool errSerious)
{
  std::string strWord(strWordC);
  std::transform(strWord.begin(), strWord.end(), strWord.begin(), (int(*)(int)) std::tolower);
  std::string s1 = getWord(true,errSerious);
  if (s1!=strWord) {
    ostringstream ssV_2;
    ssV_2 << "Expected '" << strWord << "' and not '" << s1 << "'.";
    FlxError(errSerious,"ReadStream::getWord_11", ssV_2.str(), getCurrentPos());
  }
  return s1;
}

const std::string ReadStream::getWord (const bool lowercase, const bool errSerious, const bool numbeg, const bool dblpnt) {
  ostringstream ssV;

  // check if the first character is of type STRING
  {
    const ReadStream::InpType gnt = getNextType();
    if ( gnt!=STRING && (!numbeg || gnt!=NUMERAL)) {
      ostringstream ssV_2;
      char c = getChar(false);
      ssV_2 << "Character '" << c << "' (" << int(c) << ") not expected. (Only 'A'-'Z', 'a'-'z'";
      if (numbeg) {
        ssV_2 << ", '0'-'9'";
      }
      ssV_2 << " and '_' is valid.)";
      FlxError(errSerious,"ReadStream::getWord_1", ssV_2.str(), getCurrentPos());
    }
  }
  
  // read the word
  if (dblpnt) {
    while ( getNextType() == STRING || getNextType() == NUMERAL || whatIsNextChar()==':' ) {
      if (whatIsNextChar()==':') {
        if (whatIsNextString(2)=="::") {
          ssV << getChar(false);
          ssV << getChar(false);
        } else {
          break;
        }
      } else {
        ssV << getChar(false);
      }
    }
  } else {
    while ( getNextType() == STRING || getNextType() == NUMERAL ) {
      ssV << getChar(false);
    }
  }
  setNext();
  if (lowercase) {
    std::string s1 = ssV.str();
    std::transform(s1.begin(), s1.end(), s1.begin(), (int(*)(int)) std::tolower);
    return s1;
  } else {
    return ssV.str();
  }
}

const std::string ReadStream::getText(const bool MultiLine, const bool errSerious) {
  getChar('"',errSerious,false);
  string strText = "";
  while ( whatIsNextChar() != '"' ) {
    if ( ! MultiLine ) {
      if ( whatIsNextChar() == 10 || whatIsNextChar() == 9 ) {
        ostringstream ssV_2;
        ssV_2 << "Character (" << getChar(false) << ") not expected.";
        FlxError(errSerious,"ReadStream::getText_1", ssV_2.str(), getCurrentPos());
      } 
    }
    char c = getChar(false);
    if (c==92) {
      const char cn = getChar(false);
      switch (cn) {
        case '"':
          c = '"';
          break;
        case '\\':
          c = '\\';
          break;
        case 't':
          c = '\t';
          break;
        case 'n':        // needs to come right before default:!!!
          if (MultiLine) {
            c = '\n';
            break;
          }
        default:
        {
          ostringstream ssV_2;
          ssV_2 << "Character (" << cn << ") not allowed after '\\'.";
          FlxError(errSerious,"ReadStream::getText_2", ssV_2.str(), getCurrentPos());
        }
      }
    }
    strText += c;
  }
  getChar('"',errSerious);
  setNext();
  return strText;
}

const bool ReadStream::getBool(const bool errSerious) {
  bool b1;
  if ( nextCanBeNumber() ) {
    b1 = (fabs(get_Double(errSerious))> GlobalVar.TOL());
  } else {
    if ( getNextType() != ReadStream::STRING ) {
      ostringstream ssV;
      ssV << "Expected 'true' or 'false'.";
      FlxError(errSerious,"ReadStream::getBool_1", ssV.str(), getCurrentPos()); b1 = false;
    } else {
      string str1 = getWord(true,errSerious);
      if ( str1 == "true" ) {
        b1 = true;
      } else if ( str1 == "false" ) {
        b1 = false;
      } else {
        ostringstream ssV;
        ssV << "Expected 'true' or 'false'.";
        FlxError(errSerious,"ReadStream::getBool_2", ssV.str(), getCurrentPos()); b1 = false;
      }
    }
  }
  return b1;
}


const bool ReadStream::nextCanBeNumber() {
  if ( theStream->eof() ) return false;
  char ch = theStream->peek();
  // check if the first character is of type NUMERAL or if it is a leading sign
  if ( getType( ch ) != NUMERAL && ch != '+' && ch != '-' && ch != '.' ) return false;
  else return true;
}

time_t ReadStream::getDate(const bool errSerious)
{
  // read the date from the input stream
    tuint iyear = get_UInt<tuint>(errSerious,false);
    const bool year_first = (whatIsNextChar()=='-');
    getChar(year_first?'-':'.',errSerious,false);
    const tuint imonth = get_UInt<tuint>(errSerious,false);
    getChar(year_first?'-':'.',errSerious,false);
    tuint iday = get_UInt<tuint>(errSerious);
    if (!year_first) {
      tuint ti = iday;
      iday = iyear;
      iyear = ti;
    }
  // transform to time_t
    // perform some checks
      if (imonth<1||imonth>12) {
        ostringstream ssV;
        ssV << "Date(month): Expected a value between 1 (January) and 12 (December); and not '" << imonth << "'.";
        FlxError(errSerious,"ReadStream::getDate_1", ssV.str(), getCurrentPos());
      }
      if (iday<1||iday>31) {
        ostringstream ssV;
        ssV << "Date(day): Expected a value between 1 and 31; and not '" << iday << "'.";
        FlxError(errSerious,"ReadStream::getDate_2", ssV.str(), getCurrentPos());
      }
      if (iyear<1900) {
        ostringstream ssV;
        ssV << "Date(year): The year must be at least 1900; and not '" << iyear << "'.";
        FlxError(errSerious,"ReadStream::getDate_3", ssV.str(), getCurrentPos());
      }
      tuint year = iyear;
        year -= 1900;
      const tuint month = imonth - 1;
      const tuint day = iday;
    struct tm newd = {0};
    newd.tm_hour = 0; newd.tm_min = 0; newd.tm_sec = 0;
    newd.tm_year = year;
    newd.tm_mon = month;  
    newd.tm_mday = day;
  return mktime(&newd);
}

time_t ReadStream::getTime(struct tm date_tm, const bool errSerious)
{
  // read the time
    date_tm.tm_hour = get_UInt<tuint>(errSerious,false);
    getChar(':',errSerious,false);
    date_tm.tm_min = get_UInt<tuint>(errSerious,false);
    if (whatIsNextChar()==':') {
      getChar(':',errSerious,false);
      date_tm.tm_sec = get_UInt<tuint>(errSerious,false);
    } else {
      date_tm.tm_sec = 0;
    }
    setNext();
  // transform the time
    return mktime(&date_tm);
}

const tdouble ReadStream::get_Double(const bool errSerious)
{
  bool read_number = false;
  // read the leading sign
    bool is_negative = false;
    char ch = theStream->peek();
    if (ch=='+') {
      getChar(false);
      ch = theStream->peek();
    } else if (ch=='-') {
      is_negative = true;
      getChar(false);
      ch = theStream->peek();
    }
  // read the number before the dot
    tdouble d = ZERO;        // the current number retrieved
    while (ch>=48&&ch<=57) {
      d *= 10;
      d += tdouble(ch-48);
      read_number = true;
      getChar(false);
      ch = theStream->peek();
    }
  // handel the dot
    if (ch=='.') {
      getChar(false);
      ch = theStream->peek();
      tdouble dd = ONE;
      while (ch>=48&&ch<=57) {
        dd /= 10;
        d += tdouble(ch-48)*dd;
        read_number = true;
        getChar(false);
        ch = theStream->peek();
      }
    }
  // check for errors
    if (!read_number) {
      // check for nan
        ch = theStream->peek();
        if (ch=='n') {
          getChar(false);
          ch = theStream->peek();
          if (ch=='a') {
            getChar(false);
            ch = theStream->peek();
            if (ch=='n') {
              getChar(false);
              if (getType(theStream->peek()==NONE)) {
                d = std::numeric_limits<tdouble>::quiet_NaN();
                if (is_negative) d = -d;
                setNext();
                return d;
              } else {
                theStream->putback('n');
                theStream->putback('a');
                theStream->putback('n');
              }
            } else {
              theStream->putback('a');
              theStream->putback('n');
            }
          } else {
            theStream->putback('n');
          }
        } else if (ch=='i') {
          getChar(false);
          ch = theStream->peek();
          if (ch=='n') {
            getChar(false);
            ch = theStream->peek();
            if (ch=='f') {
              getChar(false);
              if (getType(theStream->peek()==NONE)) {
                d = std::numeric_limits<tdouble>::infinity();
                if (is_negative) d = -d;
                setNext();
                return d;
              } else {
                theStream->putback('f');
                theStream->putback('n');
                theStream->putback('i');
              }
            } else {
              theStream->putback('n');
              theStream->putback('i');
            }
          } else {
            theStream->putback('i');
          }
        }
      std::ostringstream ssV;
      ssV << "A (floating point) number is required at this point.";
      FlxError(errSerious,"ReadStream::getDouble", ssV.str(), getCurrentPos());
    }
  // handel the exponent
    if (ch=='E'||ch=='e') {
      getChar(false);
      ch = theStream->peek();
      // read the sign
        bool is_negative_E = false;
        char ch = theStream->peek();
        if (ch=='+') {
          getChar(false);
          ch = theStream->peek();
        } else if (ch=='-') {
          is_negative_E = true;
          getChar(false);
          ch = theStream->peek();
        }
      // read the number
        int ie = 0;
        while (ch>=48&&ch<=57) {
          ie *= 10;
          ie += int(ch-48);
          getChar(false);
          ch = theStream->peek();
        }
        if (is_negative_E) ie = -ie;
      // compute ...
        d *= pow(tdouble(10),ie);
    }
  // change sign?
    if (is_negative) d = -d;
  // check if given in percent
    char strcirc[] = "°"; 
    if ( ch == '%' ) {                        // check if there is a '%' after the number
      getChar(false);
      d /= 100;
    } else if ( ch == strcirc[0] ) {                // check if there is a '°' after the number
      if (whatIsNextString(2) == "°") {
        getExpr("°");
        d *= (PI / 180);
      }
    }
  // set reader to next ...
    setNext();
  return d;
}

const tdouble ReadStream::get_GermanFloat(const bool errSerious)
{
    bool read_number = false;
  // read the leading sign
    bool is_negative = false;
    char ch = theStream->peek();
    if (ch=='+') {
      getChar(false);
      ch = theStream->peek();
    } else if (ch=='-') {
      is_negative = true;
      getChar(false);
      ch = theStream->peek();
    }
  // read the number before the dot
    tdouble d = ZERO;        // the current number retrieved
    tuint c;
    bool b;
    do {
      b = false;
      c = 0;
      while (ch>=48&&ch<=57) {
        d *= 10;
        d += tdouble(ch-48);
        read_number = true;
        getChar(false);
        ch = theStream->peek();
        ++c;
      }
      if (c==3) {
        if (ch=='.') {
          b = true;
          getChar(false);
          ch = theStream->peek();
        }
      }
    } while (b);
  // handel the dot
    if (ch==',') {
      getChar(false);
      ch = theStream->peek();
      tdouble dd = ONE;
      while (ch>=48&&ch<=57) {
        dd /= 10;
        d += tdouble(ch-48)*dd;
        read_number = true;
        getChar(false);
        ch = theStream->peek();
      }
    }
  // check for errors
    if (!read_number) {
      std::ostringstream ssV;
      ssV << "A number is required at this point.";
      FlxError(errSerious,"ReadStream::get_GermanFloat", ssV.str(), getCurrentPos());
    }
  // change sign?
    if (is_negative) d = -d;
  // set reader to next ...
    setNext();
  return d;
}


FlxIstream::FlxIstream(string nameV, bool errEOFv)
:name(nameV),errEOF(errEOFv)
{

}

void FlxIstream::reachedEOF()
{
  if (errEOF) {
    ostringstream ssV;
    ssV << "No more numbers to input. Input stream '" << name << "' is empty.";
    throw FlxException("FlxIstream::reachedEOF_1", ssV.str() ); 
  } else {
    ostringstream ssV;
    ssV << "Warning: No more numbers to input. Input stream '" << name << "' is empty.";
    GlobalVar.alert.alert("FlxIstream::reachedEOF_2", ssV.str() );
  }
}

const bool FlxIstream::get_vec(flxVec& v, tuint &lastIndex, const bool suppressErr)
{
  const tuint N = v.get_N();
  for (tuint i = 0; i<N;++i) {
    if (!get_value(v[i],suppressErr)) {
      lastIndex = i;
      return false;
    }
  }
  lastIndex = N;
  return true;
}

void FlxIstream_file::copyStream(FlxIstream* rhsB, bool errSerious, const bool err_check)
{
  FlxIstream_file* rhs = dynamic_cast<FlxIstream_file*> (rhsB);
  if (rhs==NULL) {
    std::ostringstream ssV;
    ssV << "It is not possible to replace a file-input-stream with a non-file-input-stream.";
    FlxError(errSerious,"FlxIstream_file::copyStream_1", ssV.str() );
  }
  if (err_check) {
    FlxIstream_file_combine* rhs2 = dynamic_cast<FlxIstream_file_combine*> (rhs);
    if (rhs2) {
      std::ostringstream ssV;
      ssV << "It is not possible to replace a file-input-stream with a file-combine-input-stream.";
      FlxError(errSerious,"FlxIstream_file::copyStream_2a", ssV.str() );
    }
    FlxIstream_file_binary* rhs3 = dynamic_cast<FlxIstream_file_binary*> (rhs);
    if (rhs3) {
      std::ostringstream ssV;
      ssV << "It is not possible to replace a file-input-stream with a file-binary-input-stream.";
      FlxError(errSerious,"FlxIstream_file::copyStream_2b", ssV.str() );
    }
  }
  #if FLX_DEBUG
    if (rhs == this) {
      std::ostringstream ssV;
      ssV << "ERROR";
      FlxError(errSerious,"FlxIstream_file::copyStream_3", ssV.str() );
    }
    if (rhs->iReader!=NULL && rhs->iReader == this->iReader) {
      std::ostringstream ssV;
      ssV << "ERROR";
      FlxError(errSerious,"FlxIstream_file::copyStream_4", ssV.str() );
    }
  #endif
  name = rhs->name;
  if (iReader) delete iReader;
    iReader = rhs->iReader;
    rhs->iReader = NULL;
  errEOF = rhs->errEOF;
  blocksize = rhs->blocksize;
  if (SeqVec) delete SeqVec;
    SeqVec = rhs->SeqVec;
    rhs->SeqVec = NULL;
  index = rhs->index;
  lastindex = rhs->lastindex;
  Cnumb = rhs->Cnumb;
  Cvec = rhs->Cvec;
  curCol = rhs->curCol;
  curVecIdx = rhs->curVecIdx;
  if (err_check) delete rhs;
}


FlxIstream_file::FlxIstream_file(string nameV, ReadStream* iReaderV, bool errEOFv, tuint blocksizeV, const tuint Cnumb, const vector< tuint >& Cvec, const bool do_read_block)
: FlxIstream(nameV,errEOFv)
, iReader(NULL), blocksize(blocksizeV), SeqVec(NULL), index(0), lastindex(0),
Cnumb(Cnumb),Cvec(Cvec),curCol(1),curVecIdx(0)
{
  if (iReader) delete iReader;
  iReader = iReaderV;
  errEOF = errEOFv;
  blocksize = blocksizeV;
  lastindex = blocksize;
  if (SeqVec) delete SeqVec;
  SeqVec = new tVec(blocksize);
  if (do_read_block) {
    try {
      read_block();
    } catch (...) {
      FLXMSG("FlxIstream_file::resetStream_2",1);
      iReader=NULL;
      if (SeqVec) { delete SeqVec; SeqVec=NULL; }
      throw;
    }
  }
}

FlxIstream_file::~FlxIstream_file()
{
  if (iReader) delete iReader;
  if (SeqVec) delete SeqVec;
}

void FlxIstream_file::set_next()
{
  const char c = iReader->whatIsNextChar();
  if (c == ',' || c==';') {
    iReader->getChar();
  }
}

void FlxIstream_file::read_block()
{
  if (Cnumb==1 && Cvec.size()==1) {
    for (tuint i = 0; i < blocksize; ++i) {
      if (iReader->check_eof()) {
        lastindex = i;
        break;
      }
      (*SeqVec)[i] = iReader->get_Double();
      if (!(iReader->check_eof())) {
        set_next();
      }
    }
  } else {
    tuint i=0;
    tdouble d;
    while (i<blocksize) {
      if ((iReader->check_eof())) {
        lastindex = i;
        break;
      }
      d = iReader->get_Double();
      if (!(iReader->check_eof())) {
        set_next();
      }
      if (curCol==Cvec[curVecIdx]) {        
        // read number
          (*SeqVec)[i] = d;
        // update Cvec-Position
          ++curVecIdx;
          ++i;
          if (curVecIdx==Cvec.size()) curVecIdx = 0;
      }
      ++curCol;
      if (curCol>Cnumb) curCol = 1;
    }
  }
  index = 0;
}

bool FlxIstream_file::get_value(tdouble& v, bool suppressErr)
{
  if (index == lastindex) {
    if (!suppressErr) reachedEOF();
    v = ZERO;
    return false;
  } else {
    v = (*SeqVec)[index];
    ++index;
    if (index == blocksize) read_block();
    return true;
  }
}



FlxIstream_file_binary::FlxIstream_file_binary(string nameV, const string& fileName, bool errEOFv, tuint blocksizeV, const tuint Cnumb, const vector< tuint >& Cvec, const bool binaryfloat)
: FlxIstream_file(nameV, NULL, errEOFv, blocksizeV, Cnumb, Cvec, false), file( new std::ifstream(fileName.c_str(), std::ios::in|std::ios::binary|std::ios::ate)), size(0), binaryfloat(binaryfloat)
{
  try {
    size = file->tellg();
    file->seekg(0,std::ios::beg);
    read_block();
  } catch (...) {
    FLXMSG("FlxIstream_file_combine::FlxIstream_file_binary",1);
    if (SeqVec) { delete SeqVec; SeqVec=NULL; }
    if (file) { delete file; file=NULL; }
    throw;
  }
}

FlxIstream_file_binary::~FlxIstream_file_binary()
{
  if (file) {
    file->close();
    delete file;
  }
}

void FlxIstream_file_binary::read_block()
{
  if (Cnumb==1 && Cvec.size()==1) {
    tdouble* sVp = &((*SeqVec)[0]);
    for (tuint i = 0; i < blocksize; ++i) {
      if (file->good()==false || file->tellg()>=size ) {
        lastindex = i;
        break;
      }
      if (binaryfloat) {
        float f;
        file->read((char *)&f,sizeof(float));
        sVp[i] = (tdouble)f;
      } else {
        file->read((char *)(sVp+i),sizeof(tdouble));
      }
    }
  } else {
    tuint i=0;
    tdouble d;
    while (i<blocksize) {
      if (file->good()==false || file->tellg()>=size ) {
        lastindex = i;
        break;
      }
      if (binaryfloat) {
        float f;
        file->read((char *)(&f),sizeof(float));
        d = (tdouble)f;
      } else {
        file->read((char *)(&d),sizeof(tdouble));
      }
      if (curCol==Cvec[curVecIdx]) {        
        // read number
          (*SeqVec)[i] = d;
        // update Cvec-Position
          ++curVecIdx;
          ++i;
          if (curVecIdx==Cvec.size()) curVecIdx = 0;
      }
      ++curCol;
      if (curCol>Cnumb) curCol = 1;
    }
  }
  index = 0;
}

void FlxIstream_file_binary::copyStream(FlxIstream* rhsB, bool errSerious, const bool err_check)
{
  FlxIstream_file_binary* rhs = dynamic_cast<FlxIstream_file_binary*> (rhsB);
  if (rhs==NULL) {
    std::ostringstream ssV;
    ssV << "It is not possible to replace a file-binary-input-stream with a non-file-binary-input-stream.";
    FlxError(errSerious,"FlxIstream_file_binary::copyStream_1", ssV.str() );
  }
  #if FLX_DEBUG
    if (rhs == this) {
      std::ostringstream ssV;
      ssV << "ERROR";
      FlxError(errSerious,"FlxIstream_file_binary::copyStream_3", ssV.str() );
    }
  #endif
  FlxIstream_file::copyStream(rhsB, errSerious,false);
  if (file) delete file;
  file = rhs->file;
  rhs->file = NULL;
  size = rhs->size;
  binaryfloat = rhs->binaryfloat;
  delete rhs;
}

const size_t FlxIstream_file_binary::get_N_numbers()
{
  if (file==NULL || file->good()==false) return 0;
  if (binaryfloat) {
    return size/sizeof(float);
  } else {
    return size/sizeof(tdouble);
  }
}


FlxIstream_file_combine::FlxIstream_file_combine(string nameV, std::vector< ReadStream* >& iReader_vec, flxVec& weight_vec, bool errEOFv, tuint blocksizeV, const tuint Cnumb, const vector< tuint >& Cvec)
: FlxIstream_file(nameV, NULL, errEOFv, blocksizeV, Cnumb, Cvec,false), iReader_vec(iReader_vec), weight_vecP(new flxVec(weight_vec))
{
   try {
    read_block();
  } catch (...) {
    FLXMSG("FlxIstream_file_combine::FlxIstream_file_combine",1);
    if (SeqVec) { delete SeqVec; SeqVec=NULL; }
    if (weight_vecP) { delete weight_vecP; weight_vecP=NULL; }
    throw;
  }
}

FlxIstream_file_combine::~FlxIstream_file_combine()
{
  if (weight_vecP) delete weight_vecP;
  for (size_t i=0;i<iReader_vec.size();++i) {
    if (iReader_vec[i]) delete iReader_vec[i];
  }
}

void FlxIstream_file_combine::copyStream(FlxIstream* rhsB, bool errSerious, const bool err_check)
{
  FlxIstream_file_combine* rhs = dynamic_cast<FlxIstream_file_combine*> (rhsB);
  if (rhs==NULL) {
    std::ostringstream ssV;
    ssV << "It is not possible to replace a file-combine-input-stream with a non-file-combine-input-stream.";
    FlxError(errSerious,"FlxIstream_file_combine::copyStream_1", ssV.str() );
  }
  #if FLX_DEBUG
    if (rhs == this) {
      std::ostringstream ssV;
      ssV << "ERROR";
      FlxError(errSerious,"FlxIstream_file_combine::copyStream_3", ssV.str() );
    }
  #endif
  FlxIstream_file::copyStream(rhsB, errSerious,false);
  for (size_t i=0;i<iReader_vec.size();++i) {
    if (iReader_vec[i]) delete iReader_vec[i];
  }
  iReader_vec = rhs->iReader_vec;
    rhs->iReader_vec.clear();
  if (weight_vecP) delete weight_vecP;
    weight_vecP = rhs->weight_vecP;
    rhs->weight_vecP = NULL;
  delete rhs;
}

void FlxIstream_file_combine::read_block()
{
  const size_t NR = iReader_vec.size();
  // check if readers are still working
    bool is_eof = true;
    for (size_t i=0;i<NR;++i) {
      if (iReader_vec[i]->check_eof()==false) {
        is_eof = false;
      }
    }
  tuint i=0;
  while (i<blocksize) {
    if (is_eof) {
      lastindex = i;
      break;
    }
    // get next value
      is_eof = true;
      tdouble sum_d = ZERO;  // sum of numbers
      tdouble sum_w = ZERO;  // sum of weights
      for (size_t j=0;j<NR;++j) {
        ReadStream_ptr& rp = iReader_vec[j];
        if (rp->check_eof()==false) {
          const tdouble w = weight_vecP->operator[](j);
          sum_d += w * (rp->get_Double());
          sum_w += w;
          if (rp->check_eof()==false) {
            const char c = rp->whatIsNextChar();
            if (c == ',' || c==';') {
              rp->getChar();
            }
            if (rp->check_eof()==false) is_eof = false;
          }
        }
      }
    if (curCol==Cvec[curVecIdx]) {        
      // read number
        #if FLX_DEBUG
          if (sum_w==ZERO) {
            throw FlxException_Crude("FlxIstream_file_combine::read_block");
          }
        #endif
        (*SeqVec)[i] = sum_d/sum_w;
      // update Cvec-Position
        ++curVecIdx;
        ++i;
        if (curVecIdx==Cvec.size()) curVecIdx = 0;
    }
    ++curCol;
    if (curCol>Cnumb) curCol = 1;
  }
  index = 0;
}


FlxIstream_vector::FlxIstream_vector(string nameV, FlxIstream* is, bool errEOFv, const tulong Nreserve)
: FlxIstream(nameV, errEOFv), index(0), numbEl(0)
{
  iVec.reserve(Nreserve);
  if (is!=NULL) {
    tdouble d;
    while (is->get_value(d,true)) {
      iVec.push_back(d);
    }
    numbEl = iVec.size();
  }
}

void FlxIstream_vector::copyStream(FlxIstream* rhsB, bool errSerious, const bool err_check)
{
  FlxIstream_vector* rhs = dynamic_cast<FlxIstream_vector*> (rhsB);
  if (rhs==NULL) {
    std::ostringstream ssV;
    ssV << "It is not possible to replace a vector-input-stream with a non-vector-input-stream.";
    FlxError(errSerious,"FlxIstream_vector::copyStream_1", ssV.str() );
  }
  #if FLX_DEBUG
    if (rhs == this) {
      std::ostringstream ssV;
      ssV << "ERROR";
      FlxError(errSerious,"FlxIstream_vector::copyStream_2", ssV.str() );
    }
  #endif
  name = rhs->name;
  errEOF = rhs->errEOF;
  index = rhs->index;
  numbEl = rhs->numbEl;
  iVec = rhs->iVec;
  delete rhs;
}

bool FlxIstream_vector::get_value(tdouble& v, bool suppressErr)
{
  if (index == numbEl) {
    if (!suppressErr) reachedEOF();
    v = ZERO;
    index = 0;
    return false;
  } else {
    v = iVec[index];
    ++index;
    return true;
  }
}

void FlxIstream_vector::appendNumber(const tdouble& d)
{
  iVec.push_back(d);
  numbEl = iVec.size();
}

void FlxIstream_vector::clear()
{
  iVec.clear();
  numbEl = iVec.size();
  index = 0;
}

void FlxIstream_vector::sortStream()
{
  std::sort(iVec.begin(),iVec.end());
  index = 0;
}



