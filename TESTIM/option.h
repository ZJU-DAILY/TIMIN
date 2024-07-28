#ifndef _OPTION_H_
#define _OPTION_H_
#include <sstream>
#include <map>

class OptionParser{

private:
	std::map<std::string, char *> paras;
	bool isValid;
public:
	OptionParser(int argc, char ** argv);
	char * getPara(const char * str);
	bool validCheck();
};

#endif
