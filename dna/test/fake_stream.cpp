#include "fake_stream.hpp"

fake_stream::fake_stream() :
		data_(0),
		chunksize_(1),
		offset_(0)
{ }

fake_stream::fake_stream(const fake_stream& other) :
		data_(other.data_),
		chunksize_(other.chunksize_),
		offset_(other.offset_.load())
{ }

fake_stream::fake_stream(fake_stream&& other) noexcept :
		data_(std::move(other.data_)),
		chunksize_(other.chunksize_),
		offset_(other.offset_.exchange(0))
{ }

fake_stream::fake_stream(std::vector<std::byte> data, std::size_t chunksize) :
		data_(std::move(data)),
		chunksize_(chunksize),
		offset_(0)
{ }

fake_stream::fake_stream(const std::string& data, std::size_t chunksize) :
		chunksize_(chunksize),
		offset_(0)
{
    data_.reserve((data.size() / 4) + 1);
    int baseCnt = 0;
    std::byte bases = static_cast<std::byte>(0);
    for (auto it = data.begin(); it != data.end(); ++it)
    {
        if (baseCnt == 4)
        {
            data_.push_back(bases);
            baseCnt = 0;
            bases = static_cast<std::byte>(0);
        } else
        {
            bases = bases >> 2;
            bases |= static_cast<std::byte>(dna::from_char(*it)) << 6;
        }
    }
}

fake_stream& fake_stream::operator=(const fake_stream& other)
{
	chunksize_ = other.chunksize_;
	data_ = other.data_;
	offset_ = other.offset_.load();

	return *this;
}

fake_stream& fake_stream::operator=(fake_stream&& other) noexcept
{
	chunksize_ = other.chunksize_;
	data_ = std::move(other.data_);
	offset_ = other.offset_.exchange(0);

	return *this;
}

void fake_stream::seek(long offset)
{
	offset_.store(std::min(std::max(offset, 0L), static_cast<long>(data_.size())));
}

long fake_stream::size() const
{
	return data_.size();
}

dna::sequence_buffer<fake_stream::byte_view> fake_stream::read()
{
	auto offset = offset_.load(std::memory_order_consume);
	while (true)
	{
		auto len = std::min(chunksize_, data_.size() - offset_);
		if (len == 0)
			return byte_view(nullptr, 0);

		if (offset_.compare_exchange_weak(offset, offset + len, std::memory_order_release))
			return byte_view(data_.data() + offset, len);
	}
}
