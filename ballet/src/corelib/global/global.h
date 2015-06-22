#ifndef __ballet_global_h
#define __ballet_global_h

template <typename T> static inline T *GetPtrHelper(T *ptr) { return ptr; }
template <typename Wrapper> static inline typename Wrapper::pointer GetPtrHelper(const Wrapper &p) { return p.data(); }

#define BALLET_DECLARE_PRIVATE(Class) \
    inline Class##Private * d_func() { return reinterpret_cast<Class##Private *>(GetPtrHelper(d_ptr)); } \
    inline const Class##Private * d_func() const { return reinterpret_cast<const Class##Private *>(GetPtrHelper(d_ptr)); } \
    friend class Class##Private;

#define BALLET_DECLARE_PUBLIC(Class) \
    inline Class * q_func() { return static_cast<Class *>(q_ptr); } \
    inline const Class * q_func() const { return static_cast<const Class *>(q_ptr); } \
    friend class Class;

#define BALLET_D(Class) Class##Private * const d = d_func()
#define BALLET_Q(Class) Class * const q = q_func()

#endif
