#ifndef LINKEDLIST_H
#define LINKEDLIST_H

using namespace std;

template<class T> class linkedList
{
friend class CBox;
protected:
    T item;
    linkedList* next;
public:
    linkedList():next(0)
    {
    } // default constructor

    linkedList(T& item, linkedList* N=0): item(item), next(N)
    {
    } // constructor

    const T& operator()() const
    {
        return item;
    } // read item field

    const linkedList* readNext() const
    {
        return next;
    } // read next

    const linkedList& operator=(const linkedList&);

    linkedList(const linkedList& list):
        item(list()),
        next(list.next ? new linkedList(*list.next) : 0)
    {
    } // copy constructor

    ~linkedList()
    {
        delete next;
        next = 0;
    } // destructor

    linkedList& last()
    {
        return next ? next->last() : *this;
    } // last item

    int length() const
    {
        return next ? next->length() + 1 : 1;
    } // number of items

    void append(T& t)
    {
        //cout << "&t = " << &t << endl;
        last().next = new linkedList(t);
    } // append a new item

    void insertNextItem(T& t)
    {
        next = new linkedList(t, next);
    } // inserts item in the second place

    void insertFirstItem(T& t)
    {
        next = new linkedList(item, next);
        item = t;
    } // insert a new item at the beginning ("push_front")

    void dropNextItem();
    void dropFirstItem();
    const linkedList& operator+=(linkedList&);
    linkedList& order(int);
};

template<class T>
const linkedList<T>&
linkedList<T>::operator=(const linkedList<T>&L)
{
    if (this != &L)
    {
        item = L();
        if (next)
        {
            if (L.next)
            {
                *next = *L.next;
            }
            else
            {
                delete next;
                next = 0;
            }
        }
        else
            if (L.next) next = new linkedList(*L.next);
    }
    return *this;
} // assignment operator

template<class T>
void linkedList<T>::dropNextItem()
{
    if (next)
    {
        if (next->next)
        {
            linkedList<T>* keep = next;
            next = next->next;
            keep->item.~T();
        }
        else
        {
            delete next;
            next = 0;
        }
    }
    else
        cout << "error: cannot drop nonexisting next item" << endl;
} // drop the second item from the linked list

template<class T>
void linkedList<T>::dropFirstItem()
{
    if (next)
    {
        item = next->item;
        dropNextItem();
    }
    else
        cout << "error: cannot drop the only item" << endl;
} // drop the first item in the linked list

//template<class T>
//void print(const linkedList<T>& list)
//{
//    cout << "item:" << endl;
//    print(list());
//    if (list.readNext()) print(*list.readNext());
//} // print a linked list recursively

//template<class T>
//const linkedList<T>&
//linkedList<T>::operator+=(linkedList<T>& L)
//{
//    linkedList<T>* runner = this;
//    linkedList<T>* Lrunner = &L;
//    if (L.item < item)
//    {
//        insertFirstItem(L.item);
//        Lrunner = L.next;
//    }
//    for (; runner->next; runner = runner->next)
//    {
//        of (Lrunner&&(Lrunner->item == runner->item))
//        {
//            runner->item += Lrunner->item;
//            Lrunner = Lrunner->next;
//        }
//        for (; Lrunner&&(Lrunner->item < runner->next->item);
//             Lrunner = Lrunner->next)
//        {
//            runner->insertNextItem(Lrunner->item);
//            runner = runner->next;
//        }
//    }
//    if (Lrunner&&(Lrunner->item == runner->item))
//    {
//        runner->item += Lrunner->item;
//        Lrunner = Lrunner->next;
//    }
//    if (Lrunner)
//        runner->next = new linkedList<T>(*Lrunner);
//    return *this;
//} // merge two linked lists while preserving order

//template<class T>
//linkedList<T>&
//linkedList<T>::order(int length)
//{
//    if (length > 1)
//    {
//        linkedList<T>* runner = this;
//        for (int i = 0; i < length/2-1; i++)
//            runner = runner->next;
//        linkedList<T>* second = runner->next;
//        runner->next = 0;
//        order(length/2);
//        *this += second->order(length-length/2);
//    }
//    return *this;
//} // order a disordered linked list

#endif // LINKEDLIST_H
