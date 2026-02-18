// If you have multiple specializations at lines 9, 13, and 24:
namespace Catch {
    template<>
    struct StringMaker<YourType1> {
        static std::string convert(YourType1 const& value) {
            // conversion logic
        }
    };
    
    template<>
    struct StringMaker<YourType2> {
        static std::string convert(YourType2 const& value) {
            // conversion logic
        }
    };
}