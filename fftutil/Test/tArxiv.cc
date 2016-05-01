#include <iostream>
#include <fstream>

#include "gcp/program/Program.h"

#include "gcp/util/CurlUtils.h"
#include "gcp/util/DirList.h"
#include "gcp/util/Exception.h"
#include "gcp/util/String.h"
#include "gcp/util/FdSet.h"
#include "gcp/util/Port.h"
#include "gcp/util/CoProc.h"

using namespace std;
using namespace gcp::util;using namespace gcp::program;void Program::initializeUsage() {};
using namespace gcp::program;

KeyTabEntry Program::keywords[] = {
  { "label",      "",                      "s", "Label to associate with this entry"},
  { "id",         "",                      "s", "Astro-ph id to fetch"},
  { "dir",        "",                      "s", "Directory to add"},
  { "tags",       "",                      "s", "Tags to associate"},

  { "incalltags", "",                      "s", "Tags to include (all)"},
  { "incanytags", "",                      "s", "Tags to include (any)"},
  { "excanytags", "",                      "s", "Tags to exclude (any)"},

  { "dbfile",     "",                      "s", "Database file"},
  { "bibtexfile", "",                      "s", "BIbtex file to parse"},
  { "add",        "f",                     "b", "True to add an entry"},
  { "list",       "f",                     "b", "True to list the database"},
  { "maxauth",    "0",                     "i", "Max numbers of authors to list"},
  { "latex",      "t",                     "b", "True to output entries as latex"},
  { "listtag",    "",                      "s", "List only entries matching listtag"},
  { "checkref",   "f",                     "b", "If true, list entries journal references"},
  { "ital",       "f",                     "b", "If true, use ital macro \\ital{} instead of \\sl for italics"},
  { END_OF_KEYWORDS}
};

//============================================================
// Class for managing an author
//============================================================

class Author {
public:
  String first_;
  String last_;

  void write(std::ofstream& fout);

  friend std::ostream& operator<<(std::ostream& os, Author& auth);
};

std::ostream& operator<<(std::ostream& os, Author& auth) {
  os << "First: '" << auth.first_ << "' Last: '" << auth.last_ << "'";
  return os;
}

void Author::write(std::ofstream& fout)
{
  fout << "<Author><First>" << first_ << "</First><Last>" << last_ << "</Last></Author>" << std::endl;
}

//============================================================
// Class for managing an entry in the database
//============================================================

class Entry {
public:

  std::vector<std::string> labels_;
  std::map<std::string, std::string> labelMap_;
  String title_;
  String booktitle_;
  String jref_;
  String link_;

  std::vector<Author> authorList_;
  std::map<std::string, std::string> authorMap_;

  std::vector<std::string> tags_;
  std::map<std::string, std::string> tagMap_;

  String journal_;
  String volume_;
  String year_;
  String number_;
  String pages_;

  Entry();
  Entry(std::string label, std::string title, std::string jref, std::vector<Author>& authors, std::string tags);
  Entry(std::string dbString);

  static std::vector<string> parseCommaSeparatedString(std::string tags);

  void addArxivLink(std::string link);
  void addTags(std::string tags);
  void parseDbTags(String& tags);
  void writeTags(std::ofstream& fout);

  void writeAuthorList(std::ofstream& fout);
  void parseDbAuthorList(String& authors);
  void write(std::ofstream& fout);

  void addAuthors(std::vector<Author>& authors);
  void addTitle(std::string title);
  void addJournal(std::string journal);

  void addTag(std::string name);
  void addKnownTags();

  bool hasTag(std::string tag);
  bool hasAnyTags(std::string tag);
  bool hasAllTags(std::string tag);

  bool matches(std::string incalltags, std::string incanytags, std::string excanytags);

  void addLabels(std::string labels);
  void parseDbLabels(String& labels);
  void addLabel(std::string name);
  bool hasLabel(std::string label);
  void writeLabels(std::ofstream& fout);

  void addYear(std::string year);
  void addPages(std::string pages);
  void addVolume(std::string volume);
  void addNumber(std::string number);

  bool hasJournalRef();
  bool isArticle();
  bool isProceeding();

  void printLatex(unsigned nAuthMax=0, bool ital=false);

  friend std::ostream& operator<<(std::ostream& os, Entry& ent);
};

void Entry::addAuthors(std::vector<Author>& authors)
{
  for(unsigned i=0; i < authors.size(); i++) {
    authorList_.push_back(authors[i]);
    authorMap_[authors[i].last_.str()] = authors[i].first_.str();
  }

  if(authors.size() < 10) {
    addTag("significant");
    return;
  }

  if(authors.size() > 1 && (authors[0].last_.contains("Leitch") || authors[1].last_.contains("Leitch")))
    addTag("significant");

  return;
}

void Entry::printLatex(unsigned nAuthMax, bool ital)
{
  if(ital)
    std::cout << "\\gitem{} {\\ital{" << title_ << ".}}" << std::endl;
  else
    std::cout << "\\gitem{} {\\sl " << title_ << ".}" << std::endl;

  unsigned nAuth = authorList_.size();

  if(nAuthMax > 0) {
    String testStr(authorList_[0].last_);
    if(testStr.contains("Collaboration"))
      ++nAuthMax;
  }

  if(nAuthMax > 0)
    nAuth = authorList_.size() > nAuthMax ? nAuthMax : authorList_.size();

  for(unsigned i=0; i < nAuth; i++) {
    Author& auth = authorList_[i];
    String lastStr(auth.last_);

    if(lastStr.contains("Collaboration")) {
      std::cout << auth.last_ << ": ";
    } else {
      std::cout << auth.first_ << " " << auth.last_;

      if(i < nAuth-1)
	std::cout <<", ";
    }
  }

  if(nAuth < authorList_.size())
    std::cout << ", et al." << std::endl;
  else
    std::cout << "." << std::endl;

  if(isArticle() && hasJournalRef()) {
    std::cout << year_ << ", " << journal_ << "\\ " << volume_;
  } else if(isArticle() && !hasJournalRef()) {
    std::cout << year_ << ", " << link_;
  } else if(isProceeding() && hasJournalRef()) {
  if(ital)
    std::cout << year_ << ", in {\\ital{" << booktitle_ << "}}\\ " << volume_;
  else
    std::cout << year_ << ", in {\\sl " << booktitle_ << "}\\ " << volume_;
  }

  if(!pages_.isEmpty())
    std::cout << pages_;

  std::cout << "." << std::endl << std::endl;
}

void Entry::addArxivLink(std::string id)
{
  String idStr(id);
  std::ostringstream link;

  if(idStr.contains("."))
    link << "arXiv:" << id;
  else
    link << "arXiv:astro-ph/" << id;
    
  link_ = link.str();
}

std::ostream& operator<<(std::ostream& os, Entry& ent) {

  os << "Entry: " << std::endl;

  for(unsigned i=0; i < ent.labels_.size(); i++) {
    os << "  Label: " << ent.labels_[i] << std::endl;
  }
  os << "  Title: " << ent.title_ << std::endl;
  os << "  Jref: "  << ent.jref_ << std::endl;

  for(unsigned i=0; i < ent.authorList_.size(); i++) {
    os << "  Author: " << ent.authorList_[i];
  }

  for(unsigned i=0; i < ent.tags_.size(); i++) {
    os << "  Tag: " << ent.tags_[i] << std::endl;
  }

  os << std::endl;
  return os;
}

Entry::Entry() {};

Entry::Entry(std::string labels, std::string title, std::string jref, std::vector<Author>& authors, std::string tags)
{
  addLabels(labels);

  addTitle(title);

  jref_  = jref;

  addAuthors(authors);

  addTags(tags);
  addKnownTags();
}

Entry::Entry(std::string dbString)
{
  String str(dbString);

  String labels    = str.findNextInstanceOf("<Labels>", true, "</Labels>", true, true);
  parseDbLabels(labels);

  title_ = str.findNextInstanceOf("<Title>", true, "</Title>", true, true).str();
  String titleStr(title_);
  titleStr.replace("&#x27;", "'");
  title_ = titleStr.str();

  booktitle_ = str.findNextInstanceOf("<Booktitle>", true, "</Booktitle>", true, true).str();
  jref_  = str.findNextInstanceOf("<Jref>",  true, "</Jref>",  true, true).str();

  link_ = str.findNextInstanceOf("<Arxiv>",  true, "</Arxiv>",  true, true);

  journal_ = str.findNextInstanceOf("<Journal>", true, "</Journal>", true, true);
  volume_  = str.findNextInstanceOf("<Volume>",  true, "</Volume>",  true, true);
  number_  = str.findNextInstanceOf("<Number>",  true, "</Number>",  true, true);
  year_    = str.findNextInstanceOf("<Year>",    true, "</Year>",    true, true);

  String authors = str.findNextInstanceOf("<Authors>", true, "</Authors>", true, true);
  parseDbAuthorList(authors);

  String tags    = str.findNextInstanceOf("<Tags>", true, "</Tags>", true, true);
  parseDbTags(tags);

  addKnownTags();
}

void Entry::addKnownTags()
{
  //------------------------------------------------------------
  // Add tags we know about
  //------------------------------------------------------------

  String titleStr(title_);
  titleStr = titleStr.toLower();

  if(titleStr.contains("spt") || (titleStr.contains("south") && titleStr.contains("pole") && titleStr.contains("telescope")))
    addTag("spt");

  if(titleStr.contains("cbi") || (titleStr.contains("cosmic") && titleStr.contains("background") && titleStr.contains("imager")))
    addTag("cbi");

  if(titleStr.contains("quad"))
    addTag("quad");

  if(titleStr.contains("keck"))
    addTag("keck");

  if(titleStr.contains("bicep2"))
    addTag("bicep2");
  else if(titleStr.contains("bicep"))
    addTag("bicep1");

  if(titleStr.contains("carma"))
    addTag("carma");

  if(titleStr.contains("polarbear"))
    addTag("polarbear");

  if(titleStr.contains("sza") || (titleStr.contains("sunyaev") && titleStr.contains("zel'dovich") && titleStr.contains("array")))
    addTag("sza");

  if(titleStr.contains("quiet"))
    addTag("quiet");

  if(titleStr.contains("quiet"))
    addTag("quiet");

  if(titleStr.contains("dasi") || (titleStr.contains("degree") && titleStr.contains("angular") && titleStr.contains("scale") && titleStr.contains("interferometer")))
    addTag("dasi");
}

bool Entry::hasJournalRef()
{
  return (isArticle() && (!journal_.isEmpty() && !volume_.isEmpty() && !year_.isEmpty())) ||
    (isProceeding() && (!booktitle_.isEmpty() && !year_.isEmpty()));
}

bool Entry::isArticle()
{
  return hasTag("article");
}

bool Entry::isProceeding()
{
  return hasTag("proceeding");
}

void Entry::addLabel(std::string label)
{
  if(labelMap_.find(label) == labelMap_.end()) {
    labels_.push_back(label);
    labelMap_[label] = label;
  }
}

bool Entry::hasLabel(std::string label)
{
  if(labelMap_.find(label) == labelMap_.end())
    return false;
  else
    return true;
}

void Entry::addJournal(std::string journal)
{
  journal_ = journal;
}

void Entry::addTitle(std::string title)
{
  title_ = title;
  addKnownTags();
}

void Entry::addTag(std::string tag)
{
  if(tagMap_.find(tag) == tagMap_.end()) {
    tags_.push_back(tag);
    tagMap_[tag] = tag;
  }
}

bool Entry::matches(std::string incalltags, std::string incanytags, std::string excanytags)
{
  String incAllStr(incalltags);
  String incAnyStr(incanytags);
  String excAnyStr(excanytags);

  if(!incAllStr.isEmpty()) {
    return hasAllTags(incalltags) && !hasAnyTags(excanytags);
  }

  if(!incAnyStr.isEmpty())
    return hasAnyTags(incanytags) && !hasAnyTags(excanytags);

  return !hasAnyTags(excanytags);
}


bool Entry::hasAllTags(std::string tag)
{
  String tagStr(tag);
  if(tagStr.isEmpty())
    return true;

  std::vector<string> taglist = parseCommaSeparatedString(tag);

  for(unsigned iTag=0; iTag < taglist.size(); iTag++) {
    if(tagMap_.find(taglist[iTag]) == tagMap_.end()) {
      return  false;
    }
  }

  return true;
}

bool Entry::hasAnyTags(std::string tag)
{
  String tagStr(tag);
  if(tagStr.isEmpty())
    return false;

  std::vector<string> taglist = parseCommaSeparatedString(tag);

  for(unsigned iTag=0; iTag < taglist.size(); iTag++) {
    if(tagMap_.find(taglist[iTag]) != tagMap_.end()) {
      return  true;
    }
  }

  return false;
}

bool Entry::hasTag(std::string tag)
{
  std::vector<string> taglist = parseCommaSeparatedString(tag);

  for(unsigned iTag=0; iTag < taglist.size(); iTag++) {
    if(tagMap_.find(tag) == tagMap_.end())
      return  false;
  }

  return true;
}

void Entry::parseDbAuthorList(String& authors)
{
  String author;

  do {
    author = authors.findNextInstanceOf("<Author>", true, "</Author>", true, true);
    if(!author.isEmpty()) {
      Author auth;
      auth.first_ = author.findNextInstanceOf("<First>", true, "</First>", true, true).str();
      auth.last_  = author.findNextInstanceOf("<Last>", true, "</Last>", true, true).str();
      authorList_.push_back(auth);
    }

  } while(!author.isEmpty());
}

void Entry::parseDbTags(String& tags)
{
  String tag;
  do {
    tag = tags.findNextInstanceOf("<Tag>", true, "</Tag>", true, true);
    if(!tag.isEmpty()) {
      addTag(tag.str());
    }

  } while(!tag.isEmpty());
}

void Entry::parseDbLabels(String& labels)
{
  String label;
  do {
    label = labels.findNextInstanceOf("<Label>", true, "</Label>", true, true);
    if(!label.isEmpty()) {
      addLabel(label.str());
    }

  } while(!label.isEmpty());
}

void Entry::addLabels(std::string labels)
{
  std::vector<string> ret = parseCommaSeparatedString(labels);
  for(unsigned i=0; i < ret.size(); i++)
    addLabel(ret[i]);
}

void Entry::addTags(std::string tags)
{
  std::vector<string> ret = parseCommaSeparatedString(tags);
  for(unsigned i=0; i < ret.size(); i++)
    addTag(ret[i]);
}

std::vector<string> Entry::parseCommaSeparatedString(std::string tags)
{
  std::vector<string> ret;
  String tagStr(tags);
  String tok;

  do {
    tok = tagStr.findNextStringSeparatedByChars(",", true);
    if(!tok.isEmpty()) {
      ret.push_back(tok.toLower().str());
    }
  } while(!tok.isEmpty());

  return ret;
}

void Entry::write(std::ofstream& fout)
{
  fout << "<Entry>"  << std::endl;

  writeLabels(fout);

  fout << "<Title>"  << title_ << "</Title>" << std::endl;

  if(!booktitle_.isEmpty())
    fout << "<Booktitle>"  << booktitle_ << "</Booktitle>" << std::endl;

  fout << "<Jref>"   << jref_  << "</Jref>"  << std::endl;

  if(!link_.isEmpty())
    fout << "<Arxiv>"   << link_  << "</Arxiv>"  << std::endl;

  if(!journal_.isEmpty())
    fout << "<Journal>"   << journal_  << "</Journal>"  << std::endl;

  if(!volume_.isEmpty())
    fout << "<Volume>"   << volume_  << "</Volume>"  << std::endl;

  if(!number_.isEmpty())
    fout << "<Number>"   << number_  << "</Number>"  << std::endl;

  if(!year_.isEmpty())
    fout << "<Year>"   << year_  << "</Year>"  << std::endl;

  writeAuthorList(fout);
  writeTags(fout);

  fout << "</Entry>" << std::endl << std::endl;
}

void Entry::writeAuthorList(std::ofstream& fout)
{
  fout << "<Authors>" << std::endl;

  for(unsigned i=0; i < authorList_.size(); i++)
    authorList_[i].write(fout);

  fout << "</Authors>" << std::endl;
}

void Entry::writeTags(std::ofstream& fout)
{
  fout << "<Tags>" << std::endl;

  for(unsigned i=0; i < tags_.size(); i++)
    fout << "<Tag>" << tags_[i] << "</Tag>" << std::endl;

  fout << "</Tags>" << std::endl;
}

void Entry::writeLabels(std::ofstream& fout)
{
  fout << "<Labels>" << std::endl;

  for(unsigned i=0; i < labels_.size(); i++)
    fout << "<Label>" << labels_[i] << "</Label>" << std::endl;

  fout << "</Labels>" << std::endl;
}

//============================================================
// Class for managing a reference database file
//============================================================

class RefDatabase {
public:

  std::vector<Entry*> entries_;
  std::map<std::string, Entry*> entryMap_;

  ~RefDatabase();

  Entry* addEntry(std::string label, std::string title, std::string jref, std::vector<Author>& authors, std::string tags);

  bool hasEntry(std::string label);
  bool hasEntry(Entry* entry);
  void addEntry(std::string entryString);
  void addEntry(Entry* entry);
  void addToMap(std::string label, Entry* entry);

  void listRefs();

  void write(std::string file);
  void load(std::string file);

  void addEntryFromArxiv(std::string id, std::string label, std::string tags);
  void addBibtexEntry(std::string entryStr);
  void addEntriesFromBibtex(std::string file);

  void addArxivEntriesFromDir(std::string dir, std::string tags);
  std::vector<Author> parseArxivAuthorList(std::string inp);
  static std::vector<Author> parseBibtexAuthorList(std::string inp);
  std::string parseArxivTitle(std::string inp);
  std::string parseArxivJref(std::string inp);
  std::string parseArxivLink(std::string inp);
  std::string getUrl(std::string url);

  void createLatex(std::string incalltags, std::string incanytags, std::string excanytags, unsigned maxauth, bool ital);

  friend std::ostream& operator<<(std::ostream& os, RefDatabase& db);
};

std::ostream& operator<<(std::ostream& os, RefDatabase& db)
{
  for(unsigned i=0; i < db.entries_.size(); i++) {
    os << *db.entries_[i];
  }

  return os;
}

void RefDatabase::createLatex(std::string incalltags, std::string incanytags, std::string excanytags, unsigned maxauth, bool ital)
{
  for(unsigned i=0; i < entries_.size(); i++) {
    Entry* entry=entries_[i];

    if(entry->matches(incalltags, incanytags, excanytags)) {
      entry->printLatex(maxauth, ital);
    }

  }
  
}

RefDatabase::~RefDatabase()
{
  for(unsigned i=0; i < entries_.size(); i++)
    delete entries_[i];
}

void RefDatabase::listRefs()
{
  for(unsigned i=0; i < entries_.size(); i++) {

    if(entries_[i]->jref_.size() == 0) {
      COUT("Article: " << entries_[i]->labels_[0] << " has ref " << entries_[i]->jref_);
      COUT("Title:   " << entries_[i]->title_);
    }
  }
}

void RefDatabase::addEntry(Entry* entry)
{
  if(!hasEntry(entry)) {
    entries_.push_back(entry);
    for(unsigned i=0; i < entry->labels_.size(); i++)
      addToMap(entry->labels_[i], entry);
  }
}

void RefDatabase::addEntry(std::string entryString)
{
  Entry* entry = new Entry(entryString);
  entries_.push_back(entry);

  for(unsigned i=0; i < entry->labels_.size(); i++)
    addToMap(entry->labels_[i], entry);
}

bool RefDatabase::hasEntry(Entry* entry)
{
  for(unsigned i=0; i < entry->labels_.size(); i++)
    if(hasEntry(entry->labels_[i]))
      return true;
  return false;
}

bool RefDatabase::hasEntry(std::string label)
{
  if(entryMap_.find(label) == entryMap_.end())
    return false;
  else
    return true;
}

Entry* RefDatabase::addEntry(std::string label, std::string title, std::string jref, std::vector<Author>& authors, std::string tags)
{
  if(entryMap_.find(label) != entryMap_.end())
    ThrowSimpleColorError("An entry for: " << label << " already exists", "red");

  Entry* entry = new Entry(label, title, jref, authors, tags);
  entries_.push_back(entry);
  addToMap(label, entry);

  return entry;
}

void RefDatabase::addToMap(std::string label, Entry* entry)
{
  entryMap_[label] = entry;
}

void RefDatabase::write(std::string file)
{
  std::ofstream fout;
  fout.open(file.c_str(), ios::out);

  if(!fout) {
    ThrowError("Unable to open file: " << file);
  }

  for(unsigned i=0; i < entries_.size(); i++)
    entries_[i]->write(fout);

  fout.close();
}

void RefDatabase::load(std::string file)
{
  std::ifstream fin;
  fin.open(file.c_str(), ios::in);

  if(!fin) {
    ThrowError("Unable to open file: " << file);
  }
  
  String line;
  ostringstream os;
  bool reading = false;

  while(!fin.eof()) {
    line.initialize();
    getline(fin, line.str());

    if(line.contains("</Entry>")) {
      os << line.str();
      reading = false;

      addEntry(os.str());
      os.str("");
    }

    if(line.contains("<Entry>")) {
      reading = true;
    }

    if(reading) {
      os << line.str();
    }
  };

  fin.close();
}

void RefDatabase::addBibtexEntry(std::string inp)
{
  Entry* entry = new Entry();

  String entryStr(inp);

  String type = entryStr.findNextInstanceOf("@", true, "{", true, false);
  if(type.contains("ARTICLE"))
    entry->addTag("article");
  if(type.contains("PROCEEDING"))
    entry->addTag("proceeding");

  String jref = entryStr.findNextInstanceOf("{", true, ",", true, false);
  entry->jref_ = jref;
  entry->addLabel(jref.str());

  String authorList = entryStr.findNextInstanceOf("author = ", true, "title =", true, false);
  authorList.replace("{{", "{");

  std::vector<Author> authors = parseBibtexAuthorList(authorList.str());
  entry->addAuthors(authors);

  entryStr.resetToBeginning();;
  entry->addTitle(entryStr.findNextInstanceOf("title = \"{", true, "}\"", true, false).str());

  entryStr.resetToBeginning();
  entry->booktitle_ = entryStr.findNextInstanceOf("booktitle = {", true, "}", true, false).str();

  entryStr.resetToBeginning();;
  entry->addJournal(entryStr.findNextInstanceOf("journal = {", true, "}", true, false).str());

  entryStr.resetToBeginning();;
  String eprint = entryStr.findNextInstanceOf("eprint = {", true, "}", true, false);

  if(!eprint.isEmpty()) {
    entry->addLabel(eprint.str());
    entry->addArxivLink(eprint.str());
  }

  entryStr.resetToBeginning();;
  String keywords = entryStr.findNextInstanceOf("keywords = {", true, "}", true, false);
  entry->addTags(keywords.str());

  entryStr.resetToBeginning();;
  entry->year_ = entryStr.findNextInstanceOf("year = ", true, ",", true, false);

  entryStr.resetToBeginning();;
  entry->volume_ = entryStr.findNextInstanceOf("volume = ", true, ",", true, false);

  entryStr.resetToBeginning();;
  entry->number_ = entryStr.findNextInstanceOf("number = ", true, ",", true, false);

  entryStr.resetToBeginning();;
  entry->pages_ = entryStr.findNextInstanceOf("pages = {", true, "}", true, false);

  if((entry->isArticle() || entry->isProceeding()) && !entry->hasJournalRef()) {
    entry->addTag("posted");
  } else {
    entry->addTag("published");
  }

  if(!hasEntry(entry)) {
    addEntry(entry);
  } else {
    delete entry;
  }
}

void RefDatabase::addEntriesFromBibtex(std::string file)
{
  std::ifstream fin;
  fin.open(file.c_str(), ios::in);

  if(!fin) {
    ThrowError("Unable to open file: " << file);
  }
  
  String line;
  ostringstream os;
  bool reading = false;

  while(!fin.eof()) {
    line.initialize();
    getline(fin, line.str());

    if(reading && line.contains("@")) {
      reading = false;
      addBibtexEntry(os.str());
      os.str("");
    }

    if(line.contains("@")) {
      reading = true;
    }

    if(reading) {
      os << line.str();
    }
  };

  // Add the last entry too

  addBibtexEntry(os.str());

  fin.close();
}

std::string RefDatabase::getUrl(std::string url)
{
  std::ostringstream os;

  os << "/usr/bin/curl " << url;
  CoProc proc(os.str());

  FILE* stdIn  = proc.stdIn()->writeFp();
  FILE* stdOut = proc.stdOut()->readFp();

  fclose(stdIn);

  char c;
  os.str("");
  while((c = (char)fgetc(stdOut)) != EOF) {
    os << c;
  }

  return os.str();
}

std::string RefDatabase::parseArxivLink(std::string inp)
{
  String str(inp);
  String ret = str.findNextInstanceOf("arxivid\">", true, "</td>", true, true);

  ret = ret.findNextInstanceOf("\">", true, "</a>", true, true);
  ret.strip('\n');
  ret.advanceToNextNonWhitespaceChar();

  return ret.str();
}

std::string RefDatabase::parseArxivJref(std::string inp)
{
  String str(inp);
  String ret = str.findNextInstanceOf("jref\">", true, "</td>", true, true);
  ret.strip('\n');
  ret.advanceToNextNonWhitespaceChar();

  return ret.str();
}

std::string RefDatabase::parseArxivTitle(std::string inp)
{
  String str(inp);
  String ret = str.findNextInstanceOf("Title:</span>", true, "</h1>", true, true);
  ret.strip('\n');
  ret.advanceToNextNonWhitespaceChar();

  return ret.str();
}

std::vector<Author> RefDatabase::parseBibtexAuthorList(std::string inp)
{
  String str(inp);
  std::vector<Author> authors;

  String author, tok;

  bool cont = true;
  std::ostringstream os;

  Author auth;
  String first, last;
  do {

    if(str.remainder().contains("{")) {

      if(str.remainder().contains("Collaboration")) {
	last  = str.findNextInstanceOf("{", true, "}", true, false);
      } else {
	last  = str.findNextInstanceOf("{", true, "},", true, false);
      }

      if(str.remainder().contains("and")) {

	if(last.contains("Collaboration")) {
	  first = str.findNextInstanceOf("}", true, "and", true, false);
	  first.strip(' ');
	} else {
	  first = str.findNextInstanceOf("},", true, "and", true, false);
	}

      } else {
	first = str.findNextInstanceOf("},", true, "},", true, false);
	cont = false;
      }
    } else {
      cont = false;
    }

    auth.last_ = last.str();
    auth.first_ = first.str();

    authors.push_back(auth);

  } while(cont);

  return authors;
}

std::vector<Author> RefDatabase::parseArxivAuthorList(std::string inp)
{
  String str(inp);
  std::vector<Author> authors;

  String author, tok;

  bool cont = true;
  std::ostringstream os;
  do {
    author = str.findNextInstanceOf("/au:+", true, "</a>", true, false);

    if(!author.isEmpty()) {
      author = author.findNextInstanceOf("\">", true, "<", false, true);
      
      std::vector<String> toks;

      do {
	tok = author.findNextString();

	if(!tok.isEmpty()) {
	  toks.push_back(tok);
	} else {

	  
	  Author auth;
	  os.str("");

	  for(unsigned i=0; i < toks.size(); i++) {

	    if(i==toks.size()-1) {
	      auth.last_ = toks[i].str();
	    } else {
	      os << toks[i].str();
	    }

	  }

	  auth.first_ = os.str();
	  authors.push_back(auth);
	}

      } while(!tok.isEmpty());

    } else {
      cont = false;
    }

  } while(cont);

  return authors;
}

void RefDatabase::addArxivEntriesFromDir(std::string dir, std::string tags)
{
  DirList dirList(dir, true);

  std::list<DirList::DirEnt> files = dirList.getFiles(false);

  for(std::list<DirList::DirEnt>::iterator entry=files.begin(); entry != files.end(); entry++) {
    String fileStr(entry->name_);

    // Get the last part of the path as the tag

    String pathStr(entry->path_);
    String tok;
    do {
      tok=pathStr.findNextInstanceOf("/", true, "/",  true, false);
    } while(!tok.isEmpty());

    String rem = pathStr.remainder();
    rem.strip('/');

    std::string id = fileStr.findNextInstanceOf(" ", false, ".pdf", true, true).str();
    
    ostringstream os;
    os << tags << "," << rem.str();

    try {

      // Add this entry, using the last part of the path as the tag

      addEntryFromArxiv(id, id, os.str());
    } catch(Exception& err) {
      COUT(err.what());
    }

  }
}

void RefDatabase::addEntryFromArxiv(std::string id, std::string label, std::string tags)
{
  String idStr(id);
  std::ostringstream url;

  if(idStr.contains("."))
    url << "http://arxiv.org/abs/" << id;
  else
    url << "http://arxiv.org/abs/astro-ph/" << id;
    
  std::string ret = getUrl(url.str());

  std::string title = parseArxivTitle(ret);
  std::string jref  = parseArxivJref(ret);
  std::vector<Author> authors = parseArxivAuthorList(ret);

  Entry* entry = addEntry(label, title, jref, authors, tags);
  entry->addTag("posted");
  entry->link_ = parseArxivLink(ret);
}

int Program::main()
{
  std::string id   = Program::getStringParameter("id");
  unsigned maxauth = Program::getIntegerParameter("maxauth");
  std::string dir  = Program::getStringParameter("dir");
  std::string tags = Program::getStringParameter("tags");
  std::string incalltags = Program::getStringParameter("incalltags");
  std::string incanytags = Program::getStringParameter("incanytags");
  std::string excanytags = Program::getStringParameter("excanytags");
  bool add         = Program::getBooleanParameter("add");
  bool list        = Program::getBooleanParameter("list");
  bool latex       = Program::getBooleanParameter("latex");
  bool checkref    = Program::getBooleanParameter("checkref");
  bool ital        = Program::getBooleanParameter("ital");
  std::string file = Program::getStringParameter("dbfile");
  std::string bibtexfile = Program::getStringParameter("bibtexfile");
  std::string listtag = Program::getStringParameter("listtag");

  std::vector<string> taglist = Entry::parseCommaSeparatedString("dasi,quad");

#if 0
  for(unsigned i=-0; i < taglist.size(); i++) {
    COUT("Found " << taglist[i]);
  }

  return 0;
#endif

  RefDatabase db;
  db.load(file);

  try {

    //============================================================
    // Adding to the database
    //============================================================

    if(add) {

      if(Program::hasValue("id")) {
	std::string label;
	
	if(Program::hasValue("label")) {
	  label = Program::getStringParameter("label");
	} else {
	  label = id;
	}
	
	db.addEntryFromArxiv(id, label, tags);

      } else if(Program::hasValue("dir")) {

	COUT("Adding dir entries");
	db.addArxivEntriesFromDir(dir, tags);

      } else if(Program::hasValue("bibtexfile")) {

	COUT("Adding entries from bibtex");
	db.addEntriesFromBibtex(bibtexfile);

      } else {
	ThrowSimpleColorError("I don't know what to add", "red");
      }
    }

    //============================================================    
    // Listing
    //============================================================    

    if(list) {
      if(latex) {
	db.createLatex(incalltags, incanytags, excanytags, maxauth, ital);
      } else {
	COUT(db);
      }

    }

    if(checkref) {
      db.listRefs();
    }

  } catch(Exception& err) {
    COUT(err.what());
  }

  if(add) {
	COUT("Adding entries from bibtex writing file");
    db.write(file);
  }

  return 0;
}
