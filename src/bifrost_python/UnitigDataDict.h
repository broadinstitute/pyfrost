#include <pybind11/pybind11.h>
#include <ColoredCDBG.hpp>

namespace py = pybind11;


#ifndef PYFROST_UNITIGDATADICT_H
#define PYFROST_UNITIGDATADICT_H

namespace pyfrost {

/**
 * Support for setting unitig data as a python dict
 */
class UnitigDataDict : public CCDBG_Data_t<UnitigDataDict> {
public:
    UnitigDataDict() { }

    UnitigDataDict(UnitigDataDict const& o) : data(o.data) { }
    UnitigDataDict(UnitigDataDict&& o) noexcept : data(std::move(o.data)) { }

    UnitigDataDict& operator=(UnitigDataDict const& o) {
        data = o.data;
        return *this;
    }

    void clear(UnitigColorMap<UnitigDataDict> const& um) {
        data.clear();
    }

    void concat(UnitigColorMap<UnitigDataDict> const& um_dest, UnitigColorMap<UnitigDataDict> const& um_src) {
        // Two unitigs get concatenated. Currently we don't merge any data because it's a new unitig.
        // Maybe in the future try to intelligently merge data? Bit hard because the dict can contain all kinds
        // of data
        data.clear();
    }

    void extract(UnitigColors const* uc_dest, UnitigColorMap<UnitigDataDict> const& um_src, bool last_extraction) {
        // Again, just clear the dictionary. Hard to handle user data in the dict that can be all kinds of types.
        data.clear();
    }

    std::string serialize(UnitigColorMap<UnitigDataDict> const& um) {
        // TODO: convert entries in the dict with two letter keys to GFA tags
        return std::string();
    }

    py::dict& getDict() {
        return data;
    }

private:
    py::dict data;

};

using PyfrostColoredUMap = UnitigColorMap<UnitigDataDict>;
using const_PyfrostColoredUMap = const_UnitigColorMap<UnitigDataDict>;

}

#endif //PYFROST_UNITIGDATADICT_H
