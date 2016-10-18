#include <iostream>
#include <stddef.h>
#include <assert.h>
#include <ipp.h>

using std::cout;
using std::endl;

class SMA
{

public:
    SMA(int period) :
            period(period), window(new Ipp32f[period]), head(NULL), tail(NULL),
            total(0)
    {
        assert(period >= 1);
    }

    ~SMA()
    {
        delete[] window;
    }

    // Adds a value to the average, pushing one out if nescessary
    void add(Ipp32f val)
    {
        // Special case: Initialization
        if (head == NULL)
        {
            head = window;
            *head = val;
            tail = head;
            inc(tail);
            total = val;
            return;
        }

        // Were we already full?
        if (head == tail)
        {
            // Fix total-cache
            total -= *head;
            // Make room
            inc(head);
        }

        // Write the value in the next spot.
        *tail = val;
        inc(tail);

        // Update our total-cache
        total += val;
    }

    // Returns the average of the last P elements added to this SMA.
    // If no elements have been added yet, returns 0.0
    Ipp32f avg(Ipp32f numer) const
    {
        ptrdiff_t size = this->size();
        if (size == 0)
        {
            return 0; // No entries => 0 average
        }
        return total * numer;
    }

    void movAvg(Ipp32f *data, int len, Ipp32f numCoefficient, Ipp32f *out)
    {

        int counter = 0;
        for (Ipp32f *itr = data; itr < (data + len); itr++)
        {
            this->add(*itr);
            //cout << "Added " << *itr << " avg: " << this->avg(numCoefficient) << endl;
            out[counter] = this->avg(numCoefficient);
            counter++;
        }

    }

private:
    int period;
    Ipp32f *window; // Holds the values to calculate the average of.

    // Logically, head is before tail
    Ipp32f *head; // Points at the oldest element we've stored.
    Ipp32f *tail; // Points at the newest element we've stored.

    Ipp32f total; // Cache the total so we don't sum everything each time.

    // Bumps the given pointer up by one.
    // Wraps to the start of the array if needed.
    void inc(Ipp32f *&p)
    {
        if (++p >= window + period)
        {
            p = window;
        }
    }

    // Returns how many numbers we have stored.
    ptrdiff_t size() const
    {
        if (head == NULL)
            return 0;
        if (head == tail)
            return period;
        return (period + tail - head) % period;
    }
};
