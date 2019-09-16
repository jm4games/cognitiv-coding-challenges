#pragma once

#include "fake_person.hpp"

class fake_person_factory
{
public:
    static fake_person new_person_with_dup_chromos();
};

class data {
public:
    static std::vector<std::byte> fake();
};
