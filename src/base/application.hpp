#pragma once

#include "logging.hpp"
#include <map>
#include <fstream>
#include <filesystem>

namespace sfem
{
    class Application
    {
    public:
        std::string name() const;

        void set_log_level(LogLevel level);

        void set_option(const std::string &name, const std::string &value);

        void log_message(const std::string &msg, LogLevel level) const;

        /// @brief Get the application instance
        static Application &instance(int argc = 0, char *argv[] = nullptr,
                                     const std::string &name = {},
                                     const std::filesystem::path &log_filename = {});

    private:
        // Constructor
        Application(int argc, char *argv[],
                    const std::string &name,
                    const std::filesystem::path &log_filename);

        // Destructor
        ~Application();

        // Disable copying
        Application(const Application &) = delete;
        Application &operator=(const Application &) = delete;

        /// @brief The application's name
        std::string name_;

        /// @brief Log level
        LogLevel log_level_;

        /// @brief Log output file
        mutable std::ofstream log_file_;

        /// @brief Options dictionary
        std::map<std::string, std::string> options_;
    };
}