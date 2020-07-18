#include "UnitigDataProxy.h"
#include "UnitigColors.h"

namespace pyfrost {

void define_UnitigDataProxy(py::module& m) {

    auto handle = py::class_<UnitigDataProxy>(m, "UnitigDataProxy")
        // Dict like interface
        .def("__getitem__", &UnitigDataProxy::getData, py::is_operator())
        .def("__setitem__", &UnitigDataProxy::setData, py::is_operator())
        .def("__delitem__", &UnitigDataProxy::delData, py::is_operator())
        .def("__contains__", &UnitigDataProxy::contains, py::is_operator())
        .def("__len__", &UnitigDataProxy::size)
        .def("__iter__", [] (UnitigDataProxy const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__str__", &UnitigDataProxy::mappedSequence)

        .def("twin", [] (UnitigDataProxy const& self) {
            return self.unitigTail().twin();
        })
        .def("full_node", [] (UnitigDataProxy const& self) {
            return UnitigDataProxy(self.mappingToFullUnitig());
        });

    // Hack to let our UnitigDataProxy inherit from the collections.abc.MutableMapping Python mixin
    auto MutableMapping = py::module::import("collections.abc").attr("MutableMapping");
    handle.attr("__bases__") = py::make_tuple(MutableMapping).attr("__add__")(handle.attr("__bases__"));
}

const std::unordered_map<std::string, UnitigMetaKeys> UnitigDataProxy::hardcoded_keys = {
    {"length",          UnitigMetaKeys::LENGTH},
    {"pos",             UnitigMetaKeys::POS},
    {"strand",          UnitigMetaKeys::STRAND},
    {"head",            UnitigMetaKeys::HEAD},
    {"tail",            UnitigMetaKeys::TAIL},
    {"mapped_sequence", UnitigMetaKeys::MAPPED_SEQUENCE},
    {"unitig_sequence", UnitigMetaKeys::UNITIG_SEQUENCE},
    {"unitig_length",   UnitigMetaKeys::UNITIG_LENGTH},
    {"is_full_mapping", UnitigMetaKeys::IS_FULL_MAPPING},
    {"colors",          UnitigMetaKeys::COLORS},
};

bool UnitigDataProxy::contains(const std::string &key) const {
    if(getMetaKey(key) == UnitigMetaKeys::NONE) {
        return getDataDict().contains(key.c_str());
    }

    return true;
}

py::object UnitigDataProxy::getData(std::string const& key) const {
    switch(getMetaKey(key)) {
        case UnitigMetaKeys::LENGTH:
            return py::cast(unitig.len);

        case UnitigMetaKeys::POS:
            return py::cast((unitig.strand || unitig.isFullMapping())
                ? unitig.dist
                : (unitig.size - unitig.getGraph()->getK() - unitig.dist));

        case UnitigMetaKeys::STRAND:
            return py::cast(unitig.strand ? Strand::FORWARD : Strand::REVERSE);

        case UnitigMetaKeys::HEAD:
            return py::cast(unitigHead());

        case UnitigMetaKeys::TAIL:
            return py::cast(unitigTail());

        case UnitigMetaKeys::MAPPED_SEQUENCE:
            return py::cast(mappedSequence());

        case UnitigMetaKeys::UNITIG_SEQUENCE:
            return py::cast(unitigSequence());

        case UnitigMetaKeys::UNITIG_LENGTH:
            return py::cast(unitigLength());

        case UnitigMetaKeys::IS_FULL_MAPPING:
            return py::cast(unitig.isFullMapping());

        case UnitigMetaKeys::COLORS:
            return py::cast(UnitigColorsProxy(unitig));

        default:
        {
            auto& data_dict = getDataDict();

            if(data_dict.contains(key)) {
                return data_dict[key.c_str()];
            } else {
                throw py::key_error("Key '" + key + "' not found in unitig metadata.");
            }
        }
    }
}

void UnitigDataProxy::setData(std::string const& key, py::object const& value) {
    if(getMetaKey(key) != UnitigMetaKeys::NONE) {
        throw py::key_error("Key '" + key + "' is read only unitig metadata.");
    }

    auto& data_dict = getDataDict();
    data_dict[key.c_str()] = value;
}

void UnitigDataProxy::delData(const py::str &key) {
    if(getMetaKey(key) != UnitigMetaKeys::NONE) {
        throw py::key_error("Key '" + key.cast<std::string>() + "' is read only unitig metadata.");
    }

    auto& data_dict = getDataDict();
    PyDict_DelItem(data_dict.ptr(), key.ptr());
}

size_t UnitigDataProxy::size() const {
    return hardcoded_keys.size() + getDataDict().size();
}

std::string UnitigDataProxy::mappedSequence() const {
    return unitig.mappedSequenceToString();
}

size_t UnitigDataProxy::mappedSequenceLength() const {
    return unitig.getGraph()->getK() + unitig.len - 1;
}

std::string UnitigDataProxy::unitigSequence() const {
    return unitig.strand ? unitig.referenceUnitigToString() : reverse_complement(unitig.referenceUnitigToString());
}

size_t UnitigDataProxy::unitigLength() const {
    return unitig.size - unitig.getGraph()->getK() + 1;
}

Kmer UnitigDataProxy::unitigHead() const {
    return unitig.strand ? unitig.getUnitigHead() : unitig.getUnitigTail().twin();
}

Kmer UnitigDataProxy::unitigTail() const {
    return unitig.strand ? unitig.getUnitigTail() : unitig.getUnitigHead().twin();
}

UnitigDataKeyIterator UnitigDataProxy::begin() const {
    auto& data_dict = getDataDict();

    return {hardcoded_keys.cbegin(), hardcoded_keys.cend(),
        data_dict.begin()};
}

UnitigDataKeyIterator UnitigDataProxy::end() const {
    auto& data_dict = getDataDict();
    return {hardcoded_keys.cend(), hardcoded_keys.cend(),
            data_dict.end()};
}

UnitigMetaKeys UnitigDataProxy::getMetaKey(std::string const &key) const {
    auto key_enum = hardcoded_keys.find(key);
    if(key_enum != hardcoded_keys.end()) {
        return key_enum->second;
    }

    return UnitigMetaKeys::NONE;
}

}
