#pragma once

#include "fake_person.hpp"

class fake_person_factory
{
public:
    static fake_person new_person_with_dup_chromos();
    static fake_person new_person_with_dup_chromos_male();
};

class data {
public:
    static std::vector<std::byte> fake();
};
