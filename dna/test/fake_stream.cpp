#include "fake_stream.hpp"

fake_stream::fake_stream() :
		data_(0),
		chunksize_(1),
		offset_(0),
    len_(0)
{ }

fake_stream::fake_stream(const fake_stream& other) :
		data_(other.data_),
		chunksize_(other.chunksize_),
		offset_(other.offset_.load()),
    len_(other.len_)
{ }

fake_stream::fake_stream(fake_stream&& other) noexcept :
		data_(std::move(other.data_)),
		chunksize_(other.chunksize_),
		offset_(other.offset_.exchange(0)),
    len_(other.len_)
{ }

fake_stream::fake_stream(std::vector<std::byte> data, std::size_t chunksize) :
		data_(std::move(data)),
		chunksize_(chunksize),
		offset_(0),
    len_(data_.size() * dna::packed_size::value)

{ }

fake_stream::fake_stream(const std::string& data, std::size_t chunksize) :
		chunksize_(chunksize),
		offset_(0),
    len_(data.size())
{
    data_.reserve((data.size() / dna::packed_size::value) + 1);
    int base_cnt = 0;
    std::byte bases = static_cast<std::byte>(0);
    for (auto it = data.begin(); it != data.end(); ++it)
    {
        if (base_cnt == dna::packed_size::value)
        {
            data_.push_back(bases);
            base_cnt = 0;
            bases = static_cast<std::byte>(dna::from_char(*it));
        } else
        {
            bases = bases << 2;
            bases |= static_cast<std::byte>(dna::from_char(*it));
        }

        ++base_cnt;
    }

    if (base_cnt > 0)
    {
        if (base_cnt != dna::packed_size::value)
            bases = bases << ((dna::packed_size::value - base_cnt) * 2);
        data_.push_back(bases);
    }
}

fake_stream& fake_stream::operator=(const fake_stream& other)
{
	chunksize_ = other.chunksize_;
	data_ = other.data_;
	offset_ = other.offset_.load();
  len_ = other.len_;

	return *this;
}

fake_stream& fake_stream::operator=(fake_stream&& other) noexcept
{
	chunksize_ = other.chunksize_;
	data_ = std::move(other.data_);
	offset_ = other.offset_.exchange(0);
  len_ = other.len_;

	return *this;
}

void fake_stream::seek(long offset)
{
	offset_.store(std::min(std::max(offset, 0L), static_cast<long>(len_)));
}

long fake_stream::size() const
{
	return len_;
}

dna::sequence_buffer<fake_stream::byte_view> fake_stream::read()
{
	auto offset = offset_.load(std::memory_order_consume);
	while (true)
	{
		auto size = std::min(chunksize_, data_.size() - offset_);
		if (size == 0)
			return byte_view(nullptr, 0);

		if (offset_.compare_exchange_weak(offset, offset + size, std::memory_order_release)) {
      size_t consumed = offset * dna::packed_size::value;
      size_t len = std::min(len_ - consumed, chunksize_ * dna::packed_size::value);
			return dna::sequence_buffer(byte_view(data_.data() + offset, size), len);
    }
	}
}
