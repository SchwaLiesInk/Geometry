import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

/**
 *
 * @author Gerwyn Jones
 */
class Node<T> {

  public T data;
  public Node<T> prev;
  public Node<T> next;

  public Node(T dat) {
    data = dat;
    prev = this;
    next = this;
  }

  public Node() {
    prev = this;
    next = this;
  }

  public Node<T> Insert(Node<T> b) {
    Node<T> c = next;
    b.next = c;
    b.prev = this;
    next = b;
    c.prev = b;
    return b;
  }

  public Node<T> Remove() {
    prev.next = next;
    next.prev = prev;
    next = prev = this;
    return this;
  }

  public void Splice(Node<T> b) {
    Node<T> a = this;
    Node<T> an = a.next;
    Node<T> bn = b.next;
    a.next = bn;
    b.next = an;
    an.prev = b;
    bn.prev = a;
  }

  public boolean hasNext() {
    return next.data != null;
  }

  public T next() {
    return next.data;
  }
};

class RingList<T> {

  protected Node<T> header;
  protected Node<T> win;
  protected int length;

  public RingList() {
    header = new Node<T>();
    header.prev = header;
    header.next = header;
    win = header;
  }

  public RingList(RingList<T> copy) {
    header = new Node<T>();
    header.prev = header;
    header.next = header;
    win = header;
    T t = copy.First();
    while (t != null) {
      Append(t);
      t = copy.Next();
    }
  }

  public RingList(Collection<T> copy) {
    header = new Node<T>();
    win = header;
    header.prev = header;
    header.next = header;
    for (T t : copy) {
      Append(t);
    }
  }

  public RingList(T[] bits) {
    header = new Node<T>();
    header.prev = header;
    header.next = header;
    win = header;
    for (T i : bits) {
      Append(i);
    }
  }

  public RingList(Node<T>[] bits) {
    header = new Node<T>();
    header.prev = header;
    header.next = header;
    win = header;
    for (Node<T> i : bits) {
      Append(i.data);
    }
  }

  public Node<T>[] GetArray() {
    Node<T>[] a = new Node[length];
    T t = First();
    int i = 0;
    while (t != null) {
      a[i] = new Node(t);
      i++;
      t = Next();
    }
    return a;
  }

  public void Clear() {
    /*while (length > 0) {
     First();
     Remove();
     }*/
    header = new Node<T>();
    header.prev = header.next = header;
    win = header;
    length = 0;
  }

  public Node<T>[] Empty() {
    Node<T>[] result = new Node[length];
    int i = 0;
    while (length > 0) {
      result[i] = new Node();
      result[i].data = First();
      Remove();
      i++;
    }
    result[i] = new Node();
    result[i].data = header.data;
    header = new Node<T>();
    header.prev = header;
    header.next = header;
    win = header;
    return result;
  }

  Node<T> GetNode() {
    if (win != header) {

      return win;
    }
    return null;
  }

  public void Destroy() {
    /*while (length > 0) {
     First();
     T data = Remove();
     }*/
    header = null;
    win = null;
    length = -1;
  }

  public T get(int at) {
    return Get(at);
  }

  public T Get(int i) {
    if (i < 0) {
      return header.data;//null
    }
    T val = First();
    int j = 1;
    while (win != null) {
      if (i == j) {
        return val;
      }
      val = Next();
      j++;
    }
    return val;
  }

  public void Replace(T obj) {
    win.data = obj;
  }

  public boolean Set(T obj) {
    T val = First();
    while (win != null) {
      if (win.data.equals(obj)) {
        return true;
      }
      val = win.data;
    }
    return false;
  }

  public boolean Set(int i) {
    if (i < 0) {
      return false;
    }
    T val = First();
    int j = 0;
    while (win != null) {
      if (i == j) {
        return true;
      }
      val = Next();
      j++;
    }
    return false;
  }

  public T Insert(T val) {
    if (val == null) {
      return null;
    }

    win.Insert(new Node<T>(val));
    length++;
    return val;
  }

  public Node<T> Start() {
    win = header.next;

    return win;
  }

  public Node<T> End() {
    win = header.prev;

    return win;
  }

  public T Prepend(T val) {
    if (val == null) {
      return null;
    }

    header.Insert(new Node<T>(val));
    length++;
    return val;
  }

  public T Append(T v) {

    if (v == null) {
      return null;
    }
    header.prev.Insert(new Node<T>(v));
    length++;
    return v;
  }

  public T Add(Node<T> v) {
    if (v == null) {
      return null;
    }

    header.prev.Insert(new Node<T>(v.data));
    length++;
    return v.data;
  }

  public T Insert(T[] val) {
    for (T i : val) {

      if (i != null) {
        win.Insert(new Node<T>(i));
        length++;
      }
    }
    return win.data;
  }

  public T Prepend(T[] val) {
    for (T i : val) {

      if (i != null) {
        header.Insert(new Node<T>(i));
        length++;
      }
    }
    return win.data;
  }

  public void Append(T[] val) {
    for (T i : val) {
      if (i != null) {
        header.prev.Insert(new Node<T>(i));
        length++;
      }
    }
  }

  public RingList<T> Append(RingList<T> copy) {

    T t = copy.First();
    while (t != null) {
      Append(t);
      t = copy.Next();
    }
    return this;
  }

  public RingList<T> AppendReverse(RingList<T> copy) {

    T t = copy.Last();
    while (t != null) {
      Append(t);
      t = copy.Prev();
    }
    return this;
  }

  public RingList<T> PasteAll(RingList<T> l) {
    Node<T> a = (Node<T>) header.prev;
    a.Splice(l.header);
    length += l.length;
    l.Clear();
    return this;
  }

  public RingList<T> Append(Node<T> l) {

    Node<T> t = l;
    while (t != null) {
      Append(t.data);
      t = t.next;
    }
    return this;
  }

  public RingList<T> PrePasteAll(RingList<T> l) {
    Node<T> a = (Node<T>) header;
    a.Splice(l.header);
    length += l.length;
    l.Clear();
    return this;
  }

  public RingList<T> Prepend(RingList<T> copy) {
    T t = copy.First();
    while (t != null) {
      Prepend(t);
      t = copy.Next();
    }
    return this;
  }

  public RingList<T> PrependReverse(RingList<T> copy) {
    T t = copy.Last();
    while (t != null) {
      Prepend(t);
      t = copy.Prev();
    }
    return this;
  }

  public T Remove() {
    if (win == header) {
      return header.data;
    }
    T val = win.data;
    win = win.prev;
    win.next.Remove();
    win = win.next;
    length--;
    return val;
  }

  public T Remove(T obj) {
    if (obj == null) {
      return null;
    }

    if (win == header) {
      return header.data;
    }

    T val = First();
    while (win != null) {
      if (win.data.equals(obj)) {
        win = win.prev;
        win.next.Remove();
        val = win.data;
        win = win.next;
        break;
      }
    }
    length--;
    return val;
  }

  public T Sub(RingList<T> r, T obj) {
    return r.Remove(obj);
  }

  public T Pass(RingList<T> r, T obj) {
    T t = r.Remove(obj);
    Append(t);
    return t;
  }

  public T Paste(RingList<T> r) {
    T t = r.Remove();
    Append(t);
    return t;
  }

  public void Value(T v) {
    if (v != null) {
      if (win != header) {
        win.data = v;
      }
    }
  }

  public T Value() {
    return win.data;
  }

  public T PreValue() {
    return win.prev.data;
  }

  public Node<T> Node() {
    if (win != header) {
      return win;
    }
    return null;
  }

  public Node<T> NextNode() {
    return win.next;
  }

  public Node<T> PreNode() {
    return win.prev;
  }

  public T Get() {
    T t = win.data;
    win = win.next;
    return t;
  }

  public T Next() {
    win = win.next;
    return win.data;
  }

  public T Prev() {
    win = win.prev;
    return win.data;
  }

  public T First() {
    win = header.next;
    return win.data;
  }

  public T Last() {
    win = header.prev;
    return win.data;
  }

  public RingList<T> Add(RingList<T> r, T a) {
    r.Insert(a);
    return r;
  }

  public int Length() {
    return length;
  }

  public boolean IsFirst() {
    return win == header.next && length > 0;
  }

  public boolean IsLast() {
    return win == header.prev && length > 0;
  }

  public boolean IsHead() {
    return win == header;
  }

  public RingList<T> ArrayToList(T[] a) {
    RingList<T> s = new RingList<T>();
    for (T i : a) {
      s.Append(i);
    }
    return s;
  }

  public void CopyToArray(T[] array) {
    //array = new T[length];
    if (length > 0) {
      array[0] = First();
      for (int i = 1; i < length; i++) {
        if (this.win.next.data != null) {
          array[i] = Next();
        }
      }
    }
  }

  public void CopyToArray(Node<T>[] array) {
    //array = new T[length];
    if (length > 0) {
      array[0] = new Node();
      array[0].data = First();
      for (int i = 1; i < length; i++) {
        if (this.win.next.data != null) {
          array[i] = new Node();
          array[i].data = Next();
        }
      }
    }
  }

  public int Contains(T a) {
    //array = new T[length];
    T b = First();
    for (int i = 1; i <= length; i++) {
      if (a.equals(b)) {
        return i;
      }
      b = Next();
    }
    return -1;
  }

  public T Move(int to) {
    //array = new T[length];
    T b = First();
    for (int i = 1; i <= to; i++) {
      b = Next();
    }
    return b;
  }

  void Swap(T a, T b) {
    T s = a;
    a = b;
    b = s;
  }

  Node<T> Trim() {
    Node<T> t = win = header.prev;
    header.prev = win.prev;
    win.prev.next = header;
    t.next = null;
    t.prev = null;
    return t;
  }

  Node<T> Cut() {
    Node<T> t = win = header.next;
    header.next = win.next;
    win.next.prev = header;
    t.next = null;
    t.prev = null;
    return t;
  }
  //delegate int TCOMP<T>(T a, T b);
  /*internal static void InsertionSort(ref T[] a, TCOMP<T> compare)
   {
   int n = a.Length;
   //on line method all data can be available at any time
   for (int i = 1; i < n; i++)
   {
   T v = a[i];
   int j = i;
   while ((j > 0) && compare(v, a[i - j]) < 0)
   {
   a[i] = a[i - j];
   j--;
   }
   a[j] = v;
   }
   }
   internal static void SelectionSort(ref T[] a, TCOMP<T> compare)
   {
   int n = a.Length;
   //off line method all data must be available from start
   for (int i = 0; i < n - 1; i++)
   {
   int m = i;
   for (int j = i + 1; j < n; j++)
   {
   if (compare(a[i], a[m]) < 0)
   {
   m = j;
   }
   }
   Swap(a[i], a[m]);
   }
   }
   */
  /*internal static void SlowSort(ref T[] unsortedList, TCOMP<T> compare)
   {
   
   for (int i = 0; i < unsortedList.Length; i++) {
   for (int j = 0; j < unsortedList.Length; j++)
   {
   if (compare(unsortedList[j],unsortedList[i])==1)
   {
   T d = unsortedList[i];
   unsortedList[i] = unsortedList[j];
   unsortedList[j] = d; //break;
   }
   } 
   }
   }*/
  /*
     }*/

  public int size() {
    return length;
  }

  public boolean isEmpty() {
    return length == 0;
  }

  public boolean contains(Object o) {
    return (Contains((T) o) > 0);
  }
  public Node<T> iterator() {
    return header.next;
  }

  public T[] toArray() {
    Object[] a = new Object[length];
    T t = First();
    int i = 0;
    while (t != null) {
      a[i] = new Node();
      a[i] = t;
      i++;
      t = Next();
    }
    return (T[]) a;
  }

  public T[] toArray(Object[] a) {
    a = new Object[length];
    T t = First();
    int i = 0;
    while (t != null) {
      a[i] = new Node();
      a[i] = t;
      i++;
      t = Next();
    }
    return (T[]) a;
  }

  public boolean add(T e) {
    Append(e);
    return true;
  }

  public boolean remove(Object obj) {
    T a = this.Remove((T) obj);
    return a != null;
  }

  public boolean containsAll(Collection c) {
    return false;//c.stream().noneMatch((o) -> ((Contains((T) o) < 0)));
  }

  public boolean addAll(Collection c) {
    /*c.stream().forEach((o) -> {
     Append((T) o);
     });*/
    return true;
  }

  public boolean removeAll(Collection c) {
    boolean ok = true;
    for (Object o : c) {
      if ((Remove((T) o) == null)) {
        ok = false;
      }
    }
    return ok;
  }

  public boolean retainAll(Collection c) {
    boolean ok = false;
    Node<T> keep = new Node();
    Node<T> ret = keep;
    Node<T> n = header;
    while (n != null) {
      for (Object o : c) {
        if (o.equals(n.data)) {
          keep.data = (T) o;
          keep.next = new Node();
          keep.next.prev = keep;
          keep = keep.next;
          ok = true;
        }
      }

      n = n.next;
    }
    if (ok) {
      keep.next = null;
      header = win = ret;
    }
    return ok;
  }

  public void clear() {
    this.Clear();
  }
};

class Stack<T> {

  protected RingList<T> list;

  public Stack() {
    list = new RingList<T>();
  }

  public Stack(T[] bits) {
    list = new RingList<T>(bits);
  }

  public Stack(RingList<T> lst, boolean copy) {
    if (copy) {

      list = new RingList<T>(lst);
    } else {

      list = lst;
    }
  }

  public void ToArray(T[] array) {

    list.CopyToArray(array);
  }

  public T Get(int at) {
    return list.Get(at);
  }

  public void Push(T v) {
    list.Prepend(v);
  }

  public void Push(T[] v) {
    for (T i : v) {
      list.Prepend(i);
    }
  }

  public T Pop() {
    list.First();
    return list.Remove();
  }

  public boolean Empty() {
    return list.Length() == 0;
  }

  public int Size() {
    return list.Length();
  }

  public T Peek() {
    return list.First();
  }

  public T Top() {
    return list.First();
  }

  public T Bottom() {
    return list.Last();
  }

  public T NextToTop() {
    list.First();
    return list.Next();
  }

  public T First(Stack<T> s) {
    return s.Peek();
  }

  public T Last(Stack<T> s) {
    return s.Bottom();
  }

  public Stack<T> Add(Stack<T> s, T v) {
    s.Push(v);
    return s;
  }

  public Stack<T> Add(Stack<T> s, T[] v) {
    for (T i : v) {
      s.Push(i);
    }
    return s;
  }

  void Reverse(T[] a, int n) {
    Stack<T> s = new Stack<T>();
    for (int i = 0; i < n; i++) {
      s.Push(a[i]);
    }
    for (int i = 0; i < n; i++) {
      a[i] = s.Pop();
    }
  }
};
