#ifndef VPT_SAVABLE_H
#define VPT_SAVABLE_H

#include <string> 

namespace vpt{
  enum SavableState { SAVABLE_NEW_FILE, SAVABLE_OPENED_FILE, SAVABLE_MODIFIED_FILE}; 
  class Savable{
  protected: 
    std::string _filePath; 
    std::string _fileName; 
    SavableState _savableState; 
  public: 
    Savable(){
      _filePath=""; 
      _fileName=""; 
      _savableState = SAVABLE_NEW_FILE; 
    }

    inline void setFileName(const std::string& n) { _fileName = n; }
    inline void setFilePath(const std::string& n) { _filePath = n; }
    inline void setSavableState(SavableState ss) {_savableState = ss; }

	inline void modified(){
		if(_savableState!=SAVABLE_NEW_FILE)
			_savableState = SAVABLE_MODIFIED_FILE; 
	}

    inline std::string getFileName() { return _fileName; }
    inline std::string getFilePath() { return _filePath; }
    inline SavableState getSavableState() { return _savableState; }
  }; 
}

#endif
