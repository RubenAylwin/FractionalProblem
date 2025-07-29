#include <Msg.h>

/**
 * @brief: Read an env var.
 */
int MSG::getEnv(std::string envVar)
{
    auto var = std::getenv(envVar.c_str());

    /* If the env var does not exist, set to -1
       so that no messages will appear. */
    if (not var) {
        return -1;
    }
    try {
        return std::stoi(var);
    } catch (std::exception & ex) {
        std::cout << "Caught exception: " << ex.what() << ", when reading envvar: " << envVar << std::endl;
        return -1;
    }
}

