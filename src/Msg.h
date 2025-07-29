#ifndef MSG_MACROS
#define MSG_MACROS
#include <string>
#include <iostream>

///////////////////////////////////////
// Macro for file specific messaging //
///////////////////////////////////////

/**
 * @brief: Function that reads and env var for messaging.
 */
namespace MSG {
    int getEnv(std::string envVar);
}

/**
 * @brief: Macro to signal messaging.
 * @desc: This has to be at the top of a file. The given string is the env var
 * that describes the messaging level (should be an int).
 */
#define useMessages(envVar) static std::string _header = envVar;   \
    static int _msgLevel = MSG::getEnv(_header)

/**
 * @brief: This signals the beginning of a message.
 * @desc: Only msg(l) with l <= than the value of the env var will be printed.
 */
#define msg(level) if (level <= _msgLevel) std::cout << _header << ">> "

/**
 * @brief: Put this at the end of a message.
 */
#define endMsg std::endl

/**
 * @brief: This is to make some variables be computed only if 
 */
#define levelIf(level) if (level <= _msgLevel)

#endif
